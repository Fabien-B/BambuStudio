#include "AeroWingGenerator.hpp"

#include "../BoundingBox.hpp"
#include "../ClipperUtils.hpp"
#include "../Flow.hpp"
#include "../ExPolygon.hpp"

#include <algorithm>
#include <cmath>
#include <limits>

namespace Slic3r::AeroWing {

namespace {

bool valid_loop(const ExtrusionLoop &loop, bool allow_hole = false)
{
    if (loop.paths.empty() || (loop.loop_role() & elrAeroWing) ||
        (!allow_hole && (loop.loop_role() & elrPerimeterHole)))
        return false;
    for (const ExtrusionPath &path : loop.paths) {
        // Overhang detection may split a hole into normal and bridging wall
        // segments. Joining that contour preserves each segment's role and flow.
        const bool supported_overhang = allow_hole && path.role() == erOverhangPerimeter;
        if (path.polyline.points.size() < 2 || (is_bridge(path.role()) && !supported_overhang))
            return false;
    }
    for (size_t i = 0; i < loop.paths.size(); ++i)
        if (loop.paths[i].last_point() != loop.paths[(i + 1) % loop.paths.size()].first_point())
            return false;

    return true;
}

bool valid_tilt(double angle, double height)
{
    return std::isfinite(angle) && angle >= -89. && angle <= 89. &&
           std::isfinite(height) && height >= 0.;
}

bool valid_orientation(double angle)
{
    return std::isfinite(angle) && angle >= -180. && angle <= 180.;
}

void rotate_loop(ExtrusionLoop &loop, double radians)
{
    const double cosine = std::cos(radians), sine = std::sin(radians);
    for (ExtrusionPath &path : loop.paths)
        path.polyline.rotate(cosine, sine);
}

// Run the axis-aligned generator in a frame rotated about the object origin.
// Only commit a successful insertion, preserving the input exactly on skips
// and failures.
template<class Generate>
bool with_orientation(ExtrusionLoop &loop, double orientation_degrees, Generate generate)
{
    const double radians = orientation_degrees * PI / 180.;
    ExtrusionLoop result = loop;
    rotate_loop(result, -radians);
    if (!generate(result)) return false;
    if (result.loop_role() & elrAeroWing) {
        rotate_loop(result, radians);
        loop = std::move(result);
    }
    return true;
}

// Keep the original contour for all intersection and clearance calculations.
// Previously inserted fins are never treated as skin or trimmed by a later fin.
struct FinPosition { coord_t x; double width; };

// Find the first hole boundary across the whole two-pass fin width, not just
// its axis. This also catches a sloping boundary or a hole beside the return pass.
double first_hole_y(const Polygons &holes, double min_x, double max_x, double root_y)
{
    double first_y = std::numeric_limits<double>::infinity();
    for (const Polygon &hole : holes) {
        double min_y = std::numeric_limits<double>::infinity();
        double max_y = -std::numeric_limits<double>::infinity();
        for (size_t i = 0; i < hole.points.size(); ++i) {
            const Point &a = hole.points[i], &b = hole.points[(i + 1) % hole.points.size()];
            const double dx = double(b.x()) - double(a.x());
            double lo = 0., hi = 1.;
            if (dx == 0.) {
                if (a.x() < min_x || a.x() > max_x) continue;
            } else {
                const double t0 = (min_x - a.x()) / dx, t1 = (max_x - a.x()) / dx;
                lo = std::max(0., std::min(t0, t1));
                hi = std::min(1., std::max(t0, t1));
                if (lo > hi) continue;
            }
            const double y0 = a.y() + lo * double(b.y() - a.y());
            const double y1 = a.y() + hi * double(b.y() - a.y());
            min_y = std::min(min_y, std::min(y0, y1));
            max_y = std::max(max_y, std::max(y0, y1));
        }
        if (max_y >= root_y)
            first_y = std::min(first_y, std::max(root_y, min_y));
    }
    return first_y;
}

bool add_stiffener_at(ExtrusionLoop &loop, const Polygon &contour, const Polygons &allowed,
                      coord_t x, std::vector<bool> &is_fin, std::vector<FinPosition> &fins,
                      const Polygons &holes = {})
{
    std::vector<double> crossings;
    for (size_t i = 0; i < contour.points.size(); ++i) {
        const Point &a = contour.points[i], &b = contour.points[(i + 1) % contour.points.size()];
        if ((a.x() <= x && x < b.x()) || (b.x() <= x && x < a.x())) {
            const double t = double(x - a.x()) / double(b.x() - a.x());
            crossings.push_back(a.y() + t * double(b.y() - a.y()));
        }
    }
    std::sort(crossings.begin(), crossings.end());
    if (crossings.size() < 2 || crossings[1] <= crossings[0])
        return true;

    struct Crossing { double y; size_t path; size_t segment; };
    Crossing entry {crossings.front(), loop.paths.size(), 0};
    for (size_t p = 0; p < loop.paths.size() && entry.path == loop.paths.size(); ++p) {
        if (is_fin[p]) continue;
        const Points &points = loop.paths[p].polyline.points;
        for (size_t i = 1; i < points.size(); ++i) {
            const Point &a = points[i - 1], &b = points[i];
            // Half-open intervals count vertices once, and exclude vertical edges.
            if ((a.x() <= x && x < b.x()) || (b.x() <= x && x < a.x())) {
                const double t = double(x - a.x()) / double(b.x() - a.x());
                const double y = a.y() + t * double(b.y() - a.y());
                if (std::abs(y - entry.y) <= 2.) {
                    entry.path = p;
                    entry.segment = i - 1;
                    break;
                }
            }
        }
    }
    // The attachment may already have been replaced by a neighbouring fin.
    if (entry.path == loop.paths.size())
        return true;

    const ExtrusionPath &outer = loop.paths[entry.path];
    if (!std::isfinite(outer.width) || !std::isfinite(outer.height) ||
        outer.width <= 0. || outer.height <= 0. || outer.mm3_per_mm <= 0.)
        return false;
    const double spacing = scale_(Flow::rounded_rectangle_extrusion_spacing(outer.width, outer.height));
    const double width = scale_(outer.width);
    const Point root(x, coord_t(std::llround(entry.y)));
    const double end_y = std::min(crossings[1], first_hole_y(holes, double(x) - spacing, double(x), entry.y));
    // Keep the existing rounded tip and one extrusion width of clearance to
    // the encountered perimeter, whether it belongs to the exterior or a hole.
    const double length = end_y - entry.y - width;
    if (spacing <= 0.)
        return false;
    if (length < 3. * width)
        return true;
    // Hole contour vertices are protected paths, but are not neighbouring fins.
    for (const FinPosition &fin : fins) {
        const double separation = std::abs(double(root.x()) - double(fin.x));
        const double required = spacing + scale_(0.5 * (outer.width + fin.width));
        if (separation < required)
            return true;
    }

    const double radius = spacing / 2.;
    auto point = [&root](double dx, double dy) {
        return Point(root.x() + coord_t(std::llround(dx)), root.y() + coord_t(std::llround(dy)));
    };

    // Taper the junction and round the tip. Adjacent passes use normal flow
    // spacing instead of retracing the same segment with twice the material.
    Polyline fin;
    // The entry on the actual skin is determined below.
    fin.points = {root, point(-spacing, width)};
    constexpr int cap_segments = 12;
    for (int i = 0; i <= cap_segments; ++i) {
        const double angle = PI * (1. - double(i) / cap_segments);
        fin.points.push_back(point(radius * std::cos(angle) - radius, length - radius + radius * std::sin(angle)));
    }
    fin.points.push_back(point(0, width));
    fin.points.push_back(root);

    // Split only the chosen exterior segment and rotate to a stable attachment.
    // Do not use nearest-point splitting on the resulting self-touching loop.
    ExtrusionPaths paths;
    std::vector<bool> fin_flags;
    auto append_part = [&paths, &fin_flags](const ExtrusionPath &source, Points points) {
        points.erase(std::unique(points.begin(), points.end()), points.end());
        if (points.size() < 2)
            return;
        ExtrusionPath path = source;
        path.polyline = Polyline(std::move(points));
        paths.push_back(std::move(path));
        fin_flags.push_back(false);
    };
    const Points &points = outer.polyline.points;
    Points tail {root};
    tail.insert(tail.end(), points.begin() + entry.segment + 1, points.end());
    append_part(outer, std::move(tail));
    for (size_t i = 1; i < loop.paths.size(); ++i) {
        const size_t p = (entry.path + i) % loop.paths.size();
        paths.push_back(loop.paths[p]);
        fin_flags.push_back(is_fin[p]);
    }
    Points head(points.begin(), points.begin() + entry.segment + 1);
    head.push_back(root);
    append_part(outer, std::move(head));

    // Remove the return to root and all trailing exterior vertices closer than
    // spacing to it, including across path boundaries. Keep root fixed as the
    // distance reference, and discard paths that no longer contain a segment.
    const double spacing_squared = spacing * spacing;
    Point last_removed = root;
    ExtrusionPath connection_source = outer;
    while (!paths.empty() &&
           (paths.back().last_point() - root).cast<double>().squaredNorm() < spacing_squared) {
        if (fin_flags.back())
            return true;
        last_removed = paths.back().last_point();
        paths.back().polyline.points.pop_back();
        if (paths.back().polyline.points.size() < 2) {
            connection_source = std::move(paths.back());
            paths.pop_back();
            fin_flags.pop_back();
        }
    }
    if (paths.empty())
        return true;

    // The retained vertex and last removed vertex straddle the spacing circle.
    // Intersect their actual skin segment with that circle. A chord directly
    // toward root can leave a concave contour depending on its tessellation.
    const Point retained = paths.back().last_point();
    const Vec2d inside = (last_removed - root).cast<double>();
    const Vec2d toward_retained = (retained - last_removed).cast<double>();
    double lo = 0., hi = 1.;
    for (int i = 0; i < 48; ++i) {
        const double t = 0.5 * (lo + hi);
        if ((inside + t * toward_retained).squaredNorm() < spacing_squared)
            lo = t;
        else
            hi = t;
    }
    const Vec2d entry_offset = inside + hi * toward_retained;
    const Point skin_entry = point(entry_offset.x(), entry_offset.y());
    if (skin_entry != retained) {
        if (fin_flags.back())
            append_part(connection_source, Points {retained, skin_entry});
        else
            paths.back().polyline.points.push_back(skin_entry);
    }
    fin.points.front() = skin_entry;

    // Validate the new connection too, including on concave sections.
    // Two integer units accommodate rounding of the attachment on the edge.
    // All edits are still local: insufficient clearance leaves loop unchanged.
    if (!diff_pl(fin, allowed).empty())
        return true;

    ExtrusionPath stiffener = outer;
    stiffener.polyline = std::move(fin);
    paths.push_back(std::move(stiffener));
    fin_flags.push_back(true);
    fins.push_back({x, outer.width});
    loop.paths = std::move(paths);
    is_fin = std::move(fin_flags);
    return true;
}

} // namespace

bool add_stiffener(ExtrusionLoop &loop, double angle_degrees, double height_from_first_spiral_layer,
                   double orientation_degrees)
{
    if (!valid_orientation(orientation_degrees)) return false;
    if (orientation_degrees != 0.)
        return with_orientation(loop, orientation_degrees, [&](ExtrusionLoop &rotated) {
            return add_stiffener(rotated, angle_degrees, height_from_first_spiral_layer);
        });
    if (!valid_loop(loop) || !valid_tilt(angle_degrees, height_from_first_spiral_layer))
        return false;
    Polygon contour = loop.polygon();
    if (!contour.is_valid() || contour.area() == 0.) return false;
    contour.make_counter_clockwise();
    const BoundingBox bounds = get_extents(contour);
    const coord_t center_x = bounds.min.x() + (bounds.max.x() - bounds.min.x()) / 2;
    const double shifted_x = double(center_x) + scale_(height_from_first_spiral_layer * std::tan(angle_degrees * PI / 180.));
    if (!std::isfinite(shifted_x)) return false;
    if (shifted_x <= bounds.min.x() || shifted_x >= bounds.max.x()) return true;
    const coord_t x = coord_t(std::llround(shifted_x));
    if (x <= bounds.min.x() || x >= bounds.max.x()) return true;
    std::vector<bool> is_fin(loop.paths.size(), false);
    std::vector<FinPosition> fins;
    if (!add_stiffener_at(loop, contour, offset(contour, 2.f), x, is_fin, fins)) return false;
    if (std::find(is_fin.begin(), is_fin.end(), true) != is_fin.end())
        loop.set_loop_role(loop.loop_role() | elrAeroWing);
    return true;
}

static bool add_grids(ExtrusionLoop &loop, const Polygon &contour, const Polygons &allowed,
                     std::vector<bool> &is_fin, const std::vector<StiffenerGrid> &grids,
                     double height_from_first_spiral_layer, coord_t grid_origin_x, const Polygons &holes = {})
{
    if (grids.empty()) return false;
    const BoundingBox bounds = get_extents(contour);
    struct Placement { coord_t x; double y; size_t path; double segment; };
    std::vector<Placement> placements;
    for (const StiffenerGrid &grid : grids) {
        if (!valid_tilt(grid.angle_degrees, height_from_first_spiral_layer) ||
            !std::isfinite(grid.spacing_mm) || grid.spacing_mm < 0.1)
            return false;
        const double pitch = scale_(grid.spacing_mm);
        const double shift = scale_(height_from_first_spiral_layer * std::tan(grid.angle_degrees * PI / 180.));
        if (!std::isfinite(pitch) || !std::isfinite(shift)) return false;
        // An infinite grid covers every layer, even when its central fin has left
        // the section. Reduce its phase instead of iterating over out-of-bounds fins.
        const double phase = std::fmod(double(grid_origin_x) + shift, pitch);
        const double first_x = phase + (std::floor((double(bounds.min.x()) - phase) / pitch) + 1.) * pitch;
        for (size_t i = 0; ; ++i) {
            const double candidate_x = first_x + double(i) * pitch;
            if (candidate_x >= bounds.max.x()) break;
            const coord_t x = coord_t(std::llround(candidate_x));
            if (x <= bounds.min.x() || x >= bounds.max.x()) continue;
            Placement placement {x, 0., loop.paths.size(), 0.};
            for (size_t p = 0; p < loop.paths.size(); ++p) {
                if (is_fin[p]) continue;
                const Points &points = loop.paths[p].polyline.points;
                for (size_t j = 1; j < points.size(); ++j) {
                    const Point &a = points[j - 1], &b = points[j];
                    if ((a.x() <= x && x < b.x()) || (b.x() <= x && x < a.x())) {
                        const double t = double(x - a.x()) / double(b.x() - a.x());
                        const double y = a.y() + t * double(b.y() - a.y());
                        if (placement.path == loop.paths.size() || y < placement.y)
                            placement = {x, y, p, double(j - 1) + t};
                    }
                }
            }
            if (placement.path != loop.paths.size()) placements.push_back(placement);
        }
    }
    // Insert in perimeter traversal order, including clockwise loops. Otherwise
    // a connection on a long edge could consume a later fin's attachment.
    std::sort(placements.begin(), placements.end(), [](const Placement &a, const Placement &b) {
        return a.path != b.path ? a.path < b.path : a.segment < b.segment;
    });
    placements.erase(std::unique(placements.begin(), placements.end(), [](const Placement &a, const Placement &b) {
        return a.x == b.x;
    }), placements.end());
    std::vector<FinPosition> fins;
    for (const Placement &placement : placements) {
        if (!add_stiffener_at(loop, contour, allowed, placement.x, is_fin, fins, holes)) return false;
    }
    return true;
}

bool add_stiffeners(ExtrusionLoop &loop, double spacing_mm, double angle_degrees,
                    double height_from_first_spiral_layer, coord_t grid_origin_x, double orientation_degrees)
{
    return add_stiffeners(loop, std::vector<StiffenerGrid>{{spacing_mm, angle_degrees}},
                          height_from_first_spiral_layer, grid_origin_x, orientation_degrees);
}

bool add_stiffeners(ExtrusionLoop &loop, const std::vector<StiffenerGrid> &grids,
                    double height_from_first_spiral_layer, coord_t grid_origin_x, double orientation_degrees)
{
    if (!valid_orientation(orientation_degrees)) return false;
    if (orientation_degrees != 0.)
        return with_orientation(loop, orientation_degrees, [&](ExtrusionLoop &rotated) {
            return add_stiffeners(rotated, grids, height_from_first_spiral_layer, grid_origin_x);
        });
    if (!valid_loop(loop)) return false;
    Polygon contour = loop.polygon();
    if (!contour.is_valid() || contour.area() == 0.) return false;
    contour.make_counter_clockwise();
    ExtrusionLoop result = loop;
    std::vector<bool> is_fin(result.paths.size(), false);
    if (!add_grids(result, contour, offset(contour, 2.f), is_fin, grids,
                   height_from_first_spiral_layer, grid_origin_x)) return false;
    if (std::find(is_fin.begin(), is_fin.end(), true) != is_fin.end()) {
        result.set_loop_role(result.loop_role() | elrAeroWing);
        loop = std::move(result);
    }
    return true;
}

namespace {

struct HoleCrossing {
    size_t path;
    size_t segment;
    Point point;
};

std::vector<HoleCrossing> hole_crossings(const ExtrusionPaths &paths, coord_t x)
{
    std::vector<HoleCrossing> crossings;
    for (size_t p = 0; p < paths.size(); ++p) {
        const Points &points = paths[p].polyline.points;
        for (size_t i = 1; i < points.size(); ++i) {
            const Point &a = points[i - 1], &b = points[i];
            if ((a.x() <= x && x < b.x()) || (b.x() <= x && x < a.x())) {
                const double t = double(x - a.x()) / double(b.x() - a.x());
                crossings.push_back({p, i - 1, Point(x, coord_t(std::llround(a.y() + t * double(b.y() - a.y()))))});
            }
        }
    }
    std::sort(crossings.begin(), crossings.end(), [](const HoleCrossing &a, const HoleCrossing &b) {
        return a.point.y() < b.point.y();
    });
    return crossings;
}

// Split the exact selected contour segment; retain its flow and all vertices.
// Owners distinguish the original exterior (0), holes (1..N), and connectors (-1).
void rotate_at(const ExtrusionPaths &paths, const std::vector<int> &owners, const HoleCrossing &at,
               ExtrusionPaths &out, std::vector<int> &out_owners)
{
    auto append = [&](Points points) {
        points.erase(std::unique(points.begin(), points.end()), points.end());
        if (points.size() < 2) return;
        ExtrusionPath part = paths[at.path];
        part.polyline = Polyline(std::move(points));
        out.push_back(std::move(part));
        out_owners.push_back(owners[at.path]);
    };
    const Points &points = paths[at.path].polyline.points;
    Points tail {at.point};
    tail.insert(tail.end(), points.begin() + at.segment + 1, points.end());
    append(std::move(tail));
    for (size_t i = 1; i < paths.size(); ++i) {
        const size_t p = (at.path + i) % paths.size();
        out.push_back(paths[p]);
        out_owners.push_back(owners[p]);
    }
    Points head(points.begin(), points.begin() + at.segment + 1);
    head.push_back(at.point);
    append(std::move(head));
}

// Remove the approach to root within one lane spacing, stopping on the actual
// boundary segment. Never trim a connector or a different contour.
bool trim_attachment(ExtrusionPaths &paths, std::vector<int> &owners,
                     const Point &root, double spacing, int owner)
{
    if (paths.empty()) return false;
    const double distance_squared = spacing * spacing;
    Point last_removed = root;
    ExtrusionPath source = paths.back();
    while (!paths.empty() && (paths.back().last_point() - root).cast<double>().squaredNorm() < distance_squared) {
        if (owners.back() != owner) return false;
        last_removed = paths.back().last_point();
        paths.back().polyline.points.pop_back();
        if (paths.back().polyline.points.size() < 2) {
            source = std::move(paths.back());
            paths.pop_back();
            owners.pop_back();
        }
    }
    if (paths.empty()) return false;
    const Point retained = paths.back().last_point();
    const Vec2d inside = (last_removed - root).cast<double>();
    const Vec2d direction = (retained - last_removed).cast<double>();
    double lo = 0., hi = 1.;
    for (int i = 0; i < 48; ++i) {
        const double t = (lo + hi) * 0.5;
        if ((inside + t * direction).squaredNorm() < distance_squared) lo = t;
        else hi = t;
    }
    const Vec2d offset = inside + hi * direction;
    const Point entry(root.x() + coord_t(std::llround(offset.x())), root.y() + coord_t(std::llround(offset.y())));
    if (entry != retained) {
        if (owners.back() == owner)
            paths.back().polyline.points.push_back(entry);
        else {
            source.polyline = Polyline {retained, entry};
            paths.push_back(std::move(source));
            owners.push_back(owner);
        }
    }
    return true;
}

void reverse_paths(ExtrusionPaths &paths, std::vector<int> &owners)
{
    std::reverse(paths.begin(), paths.end());
    std::reverse(owners.begin(), owners.end());
    for (ExtrusionPath &path : paths) path.reverse();
}

} // namespace

bool connect_holes(ExtrusionLoop &outer, const std::vector<ExtrusionLoop> &holes,
                   double spacing_mm, double angle_degrees,
                   double height_from_first_spiral_layer, coord_t grid_origin_x, double orientation_degrees)
{
    return connect_holes(outer, holes, std::vector<StiffenerGrid>{{spacing_mm, angle_degrees}},
                         height_from_first_spiral_layer, grid_origin_x, orientation_degrees);
}

bool connect_holes(ExtrusionLoop &outer, const std::vector<ExtrusionLoop> &holes,
                   const std::vector<StiffenerGrid> &grids,
                   double height_from_first_spiral_layer, coord_t grid_origin_x, double orientation_degrees)
{
    if (!valid_orientation(orientation_degrees)) return false;
    if (orientation_degrees != 0.)
        return with_orientation(outer, orientation_degrees, [&](ExtrusionLoop &rotated) {
            std::vector<ExtrusionLoop> rotated_holes = holes;
            for (ExtrusionLoop &hole : rotated_holes)
                rotate_loop(hole, -orientation_degrees * PI / 180.);
            return connect_holes(rotated, rotated_holes, grids, height_from_first_spiral_layer, grid_origin_x);
        });
    if (holes.empty())
        return add_stiffeners(outer, grids, height_from_first_spiral_layer, grid_origin_x);
    if (!valid_loop(outer)) return false;

    std::vector<ExtrusionLoop> boundaries;
    boundaries.reserve(holes.size() + 1);
    boundaries.push_back(outer);
    boundaries.back().make_counter_clockwise();
    ExPolygon material(boundaries.front().polygon());
    for (const ExtrusionLoop &hole : holes) {
        if (!valid_loop(hole, true)) return false;
        Polygon polygon = hole.polygon();
        if (!polygon.is_valid() || polygon.area() == 0.) return false;
        polygon.make_clockwise();
        material.holes.push_back(std::move(polygon));
        boundaries.push_back(hole);
        boundaries.back().make_clockwise();
    }
    std::vector<const ExtrusionLoop *> originals;
    for (const ExtrusionLoop &boundary : boundaries) originals.push_back(&boundary);
    if (!material.is_valid()) return false;
    // The printable domain excludes the interiors of every hole. Two internal
    // units allow rounding at perimeter centreline attachment points.
    const Polygons allowed = offset(material, 2.f);
    ExtrusionLoop result = boundaries.front();
    std::vector<int> owners(result.paths.size(), 0);
    std::vector<bool> connected(originals.size(), false);
    connected[0] = true;
    Polygons connector_clearance;

    for (size_t remaining = holes.size(); remaining > 0; --remaining) {
        struct Join {
            size_t child = 0;
            ExtrusionPaths paths;
            std::vector<int> owners;
            Polyline outbound;
            Polyline inbound;
            double width = 0.;
            double length = 0.;
        } best;

        for (size_t child = 1; child < originals.size(); ++child) {
            if (connected[child]) continue;
            const BoundingBox bounds = get_extents(material.holes[child - 1]);
            const coord_t x = bounds.min.x() + (bounds.max.x() - bounds.min.x()) / 2;
            const auto child_crossings = hole_crossings(originals[child]->paths, x);
            if (child_crossings.size() < 2 || child_crossings.front().point.y() == child_crossings.back().point.y())
                continue;
            // Approach from +Y: attach to the upper side of the hole and look
            // for a connected boundary above it (exterior or preceding hole).
            const HoleCrossing &to = child_crossings.back();
            for (size_t parent = 0; parent < originals.size(); ++parent) {
                if (!connected[parent]) continue;
                for (const HoleCrossing &from : hole_crossings(originals[parent]->paths, x)) {
                    const double length = double(from.point.y()) - double(to.point.y());
                    if (length <= 0. || (best.child != 0 && length >= best.length)) continue;
                    const ExtrusionPath &source = originals[parent]->paths[from.path];
                    if (!std::isfinite(source.width) || !std::isfinite(source.height) ||
                        source.width <= 0. || source.height <= 0. || source.mm3_per_mm <= 0.) return false;
                    const double width = scale_(source.width);
                    // A preceding hole may provide a bridging wall segment.
                    // The connector inherits its flow, so use its round-strand spacing.
                    const double lane_spacing = scale_(is_bridge(source.role())
                        ? Flow::bridge_extrusion_spacing(source.width)
                        : Flow::rounded_rectangle_extrusion_spacing(source.width, source.height));
                    if (lane_spacing <= 0. || length <= 2. * width) continue;
                    HoleCrossing active_from = from;
                    bool found = false;
                    for (const HoleCrossing &candidate : hole_crossings(result.paths, x)) {
                        if (owners[candidate.path] == int(parent) &&
                            (candidate.point - from.point).cast<double>().norm() <= 2.) {
                            active_from = candidate;
                            found = true;
                            break;
                        }
                    }
                    if (!found) continue;
                    ExtrusionPaths parent_paths, child_paths;
                    std::vector<int> parent_owners, child_owners;
                    rotate_at(result.paths, owners, active_from, parent_paths, parent_owners);
                    if (!trim_attachment(parent_paths, parent_owners, active_from.point, lane_spacing, int(parent))) continue;
                    const ExtrusionPaths &boundary = originals[child]->paths;
                    rotate_at(boundary, std::vector<int>(boundary.size(), int(child)), to, child_paths, child_owners);
                    // Cut the beginning of the clockwise hole traversal. Arrive
                    // at the +X port and leave at the central port after the lap.
                    reverse_paths(child_paths, child_owners);
                    if (!trim_attachment(child_paths, child_owners, to.point, lane_spacing, int(child))) {
                        // The hole is too small to keep distinct entry/exit ports.
                        // Rebuild without it so it no longer blocks other joins or fins.
                        std::vector<ExtrusionLoop> remaining_holes = holes;
                        remaining_holes.erase(remaining_holes.begin() + (child - 1));
                        if (!connect_holes(outer, remaining_holes, grids,
                                           height_from_first_spiral_layer, grid_origin_x)) return false;
                        // Commit the closed hole even if no holes or grid fins remain.
                        outer.set_loop_role(outer.loop_role() | elrAeroWing);
                        return true;
                    }
                    reverse_paths(child_paths, child_owners);
                    const Polyline inbound {child_paths.back().last_point(), parent_paths.front().first_point()};
                    if (!diff_pl(inbound, allowed).empty()) continue;
                    // The exterior is CCW and holes are CW. Approaching a hole
                    // from +Y must therefore use the +X lane on the outward leg
                    // and the centre line on return. The opposite pairing crosses
                    // the contour traversal at the attachments.
                    const coord_t lane_x = x + coord_t(std::llround(lane_spacing));
                    Polyline outbound {parent_paths.back().last_point(),
                        Point(lane_x, active_from.point.y() - coord_t(std::llround(width))),
                        Point(lane_x, to.point.y() + coord_t(std::llround(width))), child_paths.front().first_point()};
                    if (!diff_pl(outbound, allowed).empty()) continue;
                    // Do not switch lanes to squeeze past another connector:
                    // preserve the winding-compatible order and reject that route.
                    if (!intersection_pl(outbound, connector_clearance).empty() ||
                        !intersection_pl(inbound, connector_clearance).empty()) continue;
                    ExtrusionPath connector = source;
                    connector.polyline = outbound;
                    parent_paths.push_back(connector);
                    parent_owners.push_back(-1);
                    for (ExtrusionPath &path : child_paths) parent_paths.push_back(std::move(path));
                    parent_owners.insert(parent_owners.end(), child_owners.begin(), child_owners.end());
                    connector.polyline = inbound;
                    parent_paths.push_back(std::move(connector));
                    parent_owners.push_back(-1);
                    best = {child, std::move(parent_paths), std::move(parent_owners),
                            std::move(outbound), inbound, source.width, length};
                }
            }
        }
        // Larger holes still require a route through the available material.
        if (best.child == 0) return false;

        const Polygons clearance = offset(Polylines {best.outbound, best.inbound}, float(scale_(best.width)));
        connector_clearance.insert(connector_clearance.end(), clearance.begin(), clearance.end());
        result.paths = std::move(best.paths);
        owners = std::move(best.owners);
        connected[best.child] = true;
    }

    // Keep the layer seam on the exterior, even when the last join attached a
    // hole to another hole and temporarily rotated the loop to that boundary.
    for (size_t p = 0; p < owners.size(); ++p) {
        if (owners[p] != 0) continue;
        ExtrusionPaths paths;
        std::vector<int> exterior_first_owners;
        rotate_at(result.paths, owners, HoleCrossing {p, 0, result.paths[p].first_point()}, paths, exterior_first_owners);
        result.paths = std::move(paths);
        owners = std::move(exterior_first_owners);
        break;
    }

    // Ordinary grid fins must avoid both holes and the dedicated connectors.
    std::vector<bool> is_fin;
    is_fin.reserve(owners.size());
    for (int owner : owners) is_fin.push_back(owner != 0);
    const Polygons grid_allowed = diff(allowed, connector_clearance);
    if (!add_grids(result, material.contour, grid_allowed, is_fin, grids,
                   height_from_first_spiral_layer, grid_origin_x, material.holes)) return false;
    result.set_loop_role(result.loop_role() | elrAeroWing);
    outer = std::move(result);
    return true;
}

} // namespace Slic3r::AeroWing
