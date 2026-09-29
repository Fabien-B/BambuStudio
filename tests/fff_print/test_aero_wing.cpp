#include <catch2/catch.hpp>

#include "libslic3r/AeroWing/AeroWingGenerator.hpp"
#include "libslic3r/ClipperUtils.hpp"
#include "libslic3r/Flow.hpp"
#include "libslic3r/ExPolygon.hpp"
#include "libslic3r/GCode/SpiralVase.hpp"
#include "libslic3r/GCode/GCodeProcessor.hpp"
#include "libslic3r/Layer.hpp"
#include "libslic3r/Line.hpp"
#include "libslic3r/Print.hpp"

#include <algorithm>
#include <map>
#include <iomanip>
#include <fstream>
#include <sstream>
#include <boost/filesystem.hpp>

using namespace Slic3r;

namespace {
DynamicPrintConfig test_config()
{
    DynamicPrintConfig config = DynamicPrintConfig::full_print_config();
    // Static defaults contain generic enum vectors without serialization maps.
    // A loaded printer preset normally supplies these; this test uses defaults.
    for (const std::string &key : config.keys()) {
        ConfigOption *option = config.option(key);
        if (auto *values = dynamic_cast<ConfigOptionEnumsGeneric *>(option))
            values->keys_map = print_config_def.get(key)->enum_keys_map;
        if (auto *values = dynamic_cast<ConfigOptionEnumsGenericNullable *>(option))
            values->keys_map = print_config_def.get(key)->enum_keys_map;
    }
    return config;
}

ExtrusionLoop make_loop(const Polygon &polygon)
{
    const Flow flow(0.45f, 0.2f, 0.4f);
    ExtrusionPath path(erExternalPerimeter, flow.mm3_per_mm(), flow.width(), flow.height());
    path.polyline = polygon.split_at_first_point();
    return ExtrusionLoop(path);
}

Polygon rectangle(double width, double height)
{
    return Polygon {{0., 0.}, {scale_(width), 0.}, {scale_(width), scale_(height)}, {0., scale_(height)}};
}

void check_closed(const ExtrusionLoop &loop)
{
    REQUIRE(!loop.paths.empty());
    for (size_t i = 0; i < loop.paths.size(); ++i) {
        REQUIRE(loop.paths[i].polyline.points.size() >= 2);
        REQUIRE(loop.paths[i].last_point() == loop.paths[(i + 1) % loop.paths.size()].first_point());
        REQUIRE(loop.paths[i].mm3_per_mm > 0.);
        for (size_t j = 1; j < loop.paths[i].polyline.points.size(); ++j)
            REQUIRE(loop.paths[i].polyline.points[j - 1] != loop.paths[i].polyline.points[j]);
    }
}

const ExtrusionPath &last_fin(const ExtrusionLoop &loop)
{
    const auto fin = std::find_if(loop.paths.rbegin(), loop.paths.rend(), [](const ExtrusionPath &path) {
        return path.polyline.points.size() == 17;
    });
    REQUIRE(fin != loop.paths.rend());
    return *fin;
}
}

TEST_CASE("AeroWing enters the fin directly from the preceding exterior vertex", "[AeroWing]")
{
    Polygon polygon = rectangle(30, 20);
    SECTION("counter-clockwise") {}
    SECTION("clockwise") { polygon.reverse(); }
    SECTION("circular section") {
        Points points;
        for (int i = 0; i < 120; ++i) {
            const double angle = 2. * PI * i / 120.;
            points.emplace_back(scale_(15. + 15. * std::cos(angle)), scale_(15. + 15. * std::sin(angle)));
        }
        polygon = Polygon(std::move(points));
    }
    ExtrusionLoop loop = make_loop(polygon);
    const double skin_length = loop.length();
    REQUIRE(AeroWing::add_stiffener(loop));
    check_closed(loop);
    REQUIRE(loop.loop_role() & elrAeroWing);
    REQUIRE(loop.length() > skin_length + scale_(15.));
    const ExtrusionPath &fin = last_fin(loop);
    const double preserved_length = loop.length() - fin.length();
    const Point &entry = fin.first_point();
    REQUIRE(entry != loop.first_point());
    REQUIRE(preserved_length + (entry - fin.last_point()).cast<double>().norm() == Approx(skin_length));
    REQUIRE(fin.polyline.points[1] != fin.last_point());
    polygon.make_counter_clockwise();
    REQUIRE(diff_pl(fin.polyline, offset(polygon, 2.f)).empty());
    REQUIRE_FALSE(AeroWing::add_stiffener(loop));
}

TEST_CASE("AeroWing tilts around Y by shifting the fin with height", "[AeroWing]")
{
    const Polygon polygon = rectangle(100, 20);
    ExtrusionLoop loop = make_loop(polygon);
    const double angle = GENERATE(-30., 0., 30.);
    const double height = GENERATE(0., 5., 10.);
    REQUIRE(AeroWing::add_stiffener(loop, angle, height));
    check_closed(loop);
    const Points &fin = last_fin(loop).polyline.points;
    REQUIRE(unscale<double>(fin.back().x()) == Approx(50. + height * std::tan(angle * PI / 180.)).margin(1e-6));
    REQUIRE(fin.back().y() == 0);
    REQUIRE(diff_pl(last_fin(loop).polyline, offset(polygon, 2.f)).empty());
    // Each layer's return pass stays parallel to Y; only its X position changes.
    const Vec2d direction = (fin[fin.size() - 3] - fin[fin.size() - 2]).cast<double>().normalized();
    REQUIRE(direction.x() == Approx(0.).margin(1e-5));
    REQUIRE(direction.y() == Approx(1.).margin(1e-5));
    if (angle == 0. || height == 0.) {
        ExtrusionLoop default_loop = make_loop(polygon);
        REQUIRE(AeroWing::add_stiffener(default_loop));
        REQUIRE(last_fin(default_loop).polyline.points == fin);
    }
}

TEST_CASE("AeroWing leaves the perimeter unchanged when the tilted axis misses the section", "[AeroWing]")
{
    ExtrusionLoop loop = make_loop(rectangle(30, 20));
    const ExtrusionLoop original = loop;
    const double angle = GENERATE(-45., 45.);
    const double height = GENERATE(15., 20., 100.);
    REQUIRE(AeroWing::add_stiffener(loop, angle, height));
    REQUIRE(loop.loop_role() == original.loop_role());
    REQUIRE(loop.paths.size() == original.paths.size());
    for (size_t i = 0; i < loop.paths.size(); ++i) {
        REQUIRE(loop.paths[i].polyline.points == original.paths[i].polyline.points);
        REQUIRE(loop.paths[i].role() == original.paths[i].role());
        REQUIRE(loop.paths[i].width == original.paths[i].width);
        REQUIRE(loop.paths[i].height == original.paths[i].height);
        REQUIRE(loop.paths[i].mm3_per_mm == original.paths[i].mm3_per_mm);
    }
    check_closed(loop);
}

TEST_CASE("AeroWing covers a section with evenly spaced fins in one closed loop", "[AeroWing]")
{
    Polygon polygon = rectangle(100, 20);
    bool segmented = false;
    SECTION("counter-clockwise") {}
    SECTION("clockwise") { polygon.reverse(); }
    SECTION("segmented perimeter") { segmented = true; }
    ExtrusionLoop loop = make_loop(polygon);
    if (segmented) {
        const ExtrusionPath prototype = loop.paths.front();
        loop.paths.clear();
        for (size_t i = 0; i < polygon.points.size(); ++i) {
            ExtrusionPath path = prototype;
            path.polyline = Polyline {polygon.points[i], polygon.points[(i + 1) % polygon.points.size()]};
            path.mm3_per_mm *= 1. + 0.1 * i;
            loop.paths.push_back(std::move(path));
        }
    }
    REQUIRE(AeroWing::add_stiffeners(loop, 20., 0., 0., scale_(50.)));
    REQUIRE(loop.loop_role() & elrAeroWing);
    check_closed(loop);
    std::vector<coord_t> roots;
    polygon.make_counter_clockwise();
    for (const ExtrusionPath &path : loop.paths) {
        // On this rectangle only the rounded fins have this many vertices.
        if (path.polyline.points.size() < 17) continue;
        roots.push_back(path.last_point().x());
        REQUIRE(diff_pl(path.polyline, offset(polygon, 2.f)).empty());
    }
    std::sort(roots.begin(), roots.end());
    REQUIRE(roots.size() == 5);
    for (size_t i = 0; i < roots.size(); ++i)
        REQUIRE(unscale<double>(roots[i]) == Approx(10. + 20. * i));
}

TEST_CASE("AeroWing keeps covering upper layers after the central tilted fin leaves", "[AeroWing]")
{
    const double angle = GENERATE(-30., 30.);
    const double height = GENERATE(0., 100., 105.);
    ExtrusionLoop loop = make_loop(rectangle(100, 20));
    REQUIRE(AeroWing::add_stiffeners(loop, 20., angle, height, scale_(50.)));
    check_closed(loop);
    std::vector<double> roots;
    for (const ExtrusionPath &path : loop.paths)
        if (path.polyline.points.size() >= 17)
            roots.push_back(unscale<double>(path.last_point().x()));
    std::sort(roots.begin(), roots.end());
    REQUIRE(roots.size() == 5);
    REQUIRE(roots.front() < 20.);
    REQUIRE(roots.back() > 80.);
    for (size_t i = 0; i < roots.size(); ++i) {
        const double grid_index = (roots[i] - 50. - height * std::tan(angle * PI / 180.)) / 20.;
        REQUIRE(grid_index == Approx(std::round(grid_index)).margin(1e-6));
        if (i > 0) REQUIRE(roots[i] - roots[i - 1] == Approx(20.).margin(1e-6));
    }
}

TEST_CASE("AeroWing combines independent grids without duplicate fins", "[AeroWing]")
{
    const double height = GENERATE(0., 5., 100.);
    const bool clockwise = GENERATE(false, true);
    const bool swap_grids = GENERATE(false, true);
    Polygon polygon = rectangle(100, 20);
    if (clockwise) polygon.reverse();
    ExtrusionLoop loop = make_loop(polygon);
    std::vector<AeroWing::StiffenerGrid> grids {{20., 45.}, {30., -45.}};
    if (swap_grids) std::reverse(grids.begin(), grids.end());
    REQUIRE(AeroWing::add_stiffeners(loop, grids, height, scale_(50.)));
    check_closed(loop);
    std::vector<double> roots;
    polygon.make_counter_clockwise();
    for (const ExtrusionPath &path : loop.paths) {
        if (path.polyline.points.size() < 17) continue;
        roots.push_back(unscale<double>(path.last_point().x()));
        REQUIRE(diff_pl(path.polyline, offset(polygon, 2.f)).empty());
    }
    std::sort(roots.begin(), roots.end());
    const std::vector<double> expected = height == 0. ? std::vector<double>{10., 20., 30., 50., 70., 80., 90.} :
                                        height == 5. ? std::vector<double>{15., 35., 45., 55., 75., 95.} :
                                                       std::vector<double>{10., 30., 40., 50., 70., 90.};
    REQUIRE(roots.size() == expected.size());
    for (size_t i = 0; i < roots.size(); ++i)
        REQUIRE(roots[i] == Approx(expected[i]).margin(2e-6));
}

TEST_CASE("AeroWing keeps clearance between nearly coincident grids", "[AeroWing]")
{
    const bool clockwise = GENERATE(false, true);
    Polygon polygon = rectangle(100, 20);
    if (clockwise) polygon.reverse();
    ExtrusionLoop loop = make_loop(polygon);
    const std::vector<AeroWing::StiffenerGrid> grids {{20., 0.}, {20., 45.}};
    REQUIRE(AeroWing::add_stiffeners(loop, grids, 0.2, scale_(50.)));
    check_closed(loop);
    Polygons occupied;
    size_t fins = 0;
    for (const ExtrusionPath &path : loop.paths) {
        if (path.polyline.points.size() < 17) continue;
        REQUIRE(intersection_pl(path.polyline, occupied).empty());
        const Polygons clearance = offset(Polylines {path.polyline}, float(scale_(path.width)));
        occupied.insert(occupied.end(), clearance.begin(), clearance.end());
        ++fins;
    }
    REQUIRE(fins == 5);
}

TEST_CASE("AeroWing rejects an invalid second grid without changing the perimeter", "[AeroWing]")
{
    ExtrusionLoop loop = make_loop(rectangle(100, 20));
    const ExtrusionLoop original = loop;
    std::vector<AeroWing::StiffenerGrid> grids {{20., 0.}, {30., 45.}};
    SECTION("invalid spacing") { grids.back().spacing_mm = 0.; }
    SECTION("invalid angle") { grids.back().angle_degrees = 90.; }
    REQUIRE_FALSE(AeroWing::add_stiffeners(loop, grids));
    REQUIRE(loop.loop_role() == original.loop_role());
    REQUIRE(loop.paths.size() == original.paths.size());
    REQUIRE(loop.paths.front().polyline.points == original.paths.front().polyline.points);
}

TEST_CASE("AeroWing retains fins across small concave skin variations", "[AeroWing]")
{
    const double bump = GENERATE(0., 0.04, 0.08);
    const bool clockwise = GENERATE(false, true);
    const bool segmented = GENERATE(false, true);
    // The middle layer's tiny inward bump used to make the straight connection
    // from the last retained vertex to root leave the skin and drop the fin.
    Polygon polygon {{0., 0.}, {scale_(14.2), 0.}, {scale_(14.8), scale_(bump)},
                     {scale_(15.), 0.}, {scale_(30.), 0.},
                     {scale_(30.), scale_(20.)}, {0., scale_(20.)}};
    if (clockwise) polygon.reverse();
    ExtrusionLoop loop = make_loop(polygon);
    if (segmented) {
        const ExtrusionPath prototype = loop.paths.front();
        loop.paths.clear();
        for (size_t i = 0; i < polygon.points.size(); ++i) {
            ExtrusionPath path = prototype;
            path.polyline = Polyline {polygon.points[i], polygon.points[(i + 1) % polygon.points.size()]};
            loop.paths.push_back(std::move(path));
        }
    }
    REQUIRE(AeroWing::add_stiffeners(loop, 10., 0., 0., scale_(15.)));
    check_closed(loop);
    polygon.make_counter_clockwise();
    std::vector<coord_t> roots;
    for (const ExtrusionPath &path : loop.paths) {
        if (path.polyline.points.size() < 17) continue;
        roots.push_back(path.last_point().x());
        REQUIRE(diff_pl(path.polyline, offset(polygon, 2.f)).empty());
        const double entry_distance = (path.first_point() - path.last_point()).cast<double>().norm();
        const double spacing = scale_(Flow::rounded_rectangle_extrusion_spacing(path.width, path.height));
        REQUIRE(entry_distance == Approx(spacing).margin(2.));
    }
    std::sort(roots.begin(), roots.end());
    REQUIRE(roots.size() == 3);
    REQUIRE(roots[1] == scale_(15.));
}

TEST_CASE("AeroWing leaves a section intact when none of the repeated fins fit", "[AeroWing]")
{
    ExtrusionLoop loop = make_loop(rectangle(100, 0.8));
    const ExtrusionLoop original = loop;
    REQUIRE(AeroWing::add_stiffeners(loop, 20.));
    REQUIRE(loop.loop_role() == original.loop_role());
    REQUIRE(loop.paths.size() == original.paths.size());
    REQUIRE(loop.paths.front().polyline.points == original.paths.front().polyline.points);
    REQUIRE(loop.total_volume() == original.total_volume());
}

TEST_CASE("AeroWing continues placing fins after one has insufficient clearance", "[AeroWing]")
{
    const Polygon polygon {{0., 0.}, {scale_(100.), 0.},
                           {scale_(100.), scale_(20.)}, {0., scale_(0.5)}};
    ExtrusionLoop loop = make_loop(polygon);
    REQUIRE(AeroWing::add_stiffeners(loop, 20., 0., 0., scale_(1.)));
    check_closed(loop);
    std::vector<coord_t> roots;
    for (const ExtrusionPath &path : loop.paths) {
        if (path.polyline.points.size() < 17) continue;
        roots.push_back(path.last_point().x());
        REQUIRE(diff_pl(path.polyline, offset(polygon, 2.f)).empty());
    }
    std::sort(roots.begin(), roots.end());
    REQUIRE(roots.size() == 4);
    for (size_t i = 0; i < roots.size(); ++i)
        REQUIRE(unscale<double>(roots[i]) == Approx(21. + 20. * i));
}

TEST_CASE("AeroWing retains a segmented skin's flow and connectivity", "[AeroWing]")
{
    const Polygon polygon = rectangle(30, 20);
    ExtrusionLoop loop = make_loop(polygon);
    const ExtrusionPath prototype = loop.paths.front();
    loop.paths.clear();
    for (size_t i = 0; i < polygon.points.size(); ++i) {
        ExtrusionPath path = prototype;
        path.polyline = Polyline {polygon.points[i], polygon.points[(i + 1) % polygon.points.size()]};
        path.mm3_per_mm *= 1. + 0.1 * i;
        loop.paths.push_back(path);
    }
    const double skin_volume = loop.total_volume();
    REQUIRE(AeroWing::add_stiffener(loop));
    check_closed(loop);
    const ExtrusionPath &fin = last_fin(loop);
    const double removed_volume = unscale<double>((fin.first_point() - fin.last_point()).cast<double>().norm()) * prototype.mm3_per_mm;
    REQUIRE(loop.total_volume() - fin.total_volume() + removed_volume == Approx(skin_volume));
}

TEST_CASE("AeroWing rejects unsupported geometry without changing its loop", "[AeroWing]")
{
    ExtrusionLoop loop = make_loop(rectangle(30, 20));
    SECTION("hole perimeter") { loop.set_loop_role(elrPerimeterHole); }
    SECTION("disconnected paths") { loop.paths.front().polyline.points.pop_back(); }
    SECTION("empty path") { loop.paths.emplace_back(erExternalPerimeter); }
    const double length = loop.length();
    REQUIRE_FALSE(AeroWing::add_stiffener(loop));
    REQUIRE(loop.length() == length);
    REQUIRE_FALSE(loop.loop_role() & elrAeroWing);
}

TEST_CASE("AeroWing preserves the original perimeter when there is insufficient room", "[AeroWing]")
{
    ExtrusionLoop loop = make_loop(rectangle(30, 0.8));
    SECTION("too short for the tip") {}
    SECTION("too narrow for the adjacent passes") { loop = make_loop(rectangle(0.3, 20)); }
    SECTION("concave section whose fin would leave the skin") {
        // A narrow neck on the centre line is too small for the round trip.
        loop = make_loop(Polygon {{0., 0.}, {scale_(30), 0.}, {scale_(30), scale_(1)},
            {scale_(15.1), scale_(1)}, {scale_(15.1), scale_(20)}, {scale_(14.9), scale_(20)},
            {scale_(14.9), scale_(1)}, {0., scale_(1)}});
    }
    const ExtrusionLoop original = loop;
    REQUIRE(AeroWing::add_stiffener(loop));
    REQUIRE(loop.loop_role() == original.loop_role());
    REQUIRE(loop.paths.size() == original.paths.size());
    for (size_t i = 0; i < loop.paths.size(); ++i) {
        REQUIRE(loop.paths[i].polyline.points == original.paths[i].polyline.points);
        REQUIRE(loop.paths[i].role() == original.paths[i].role());
        REQUIRE(loop.paths[i].width == original.paths[i].width);
        REQUIRE(loop.paths[i].height == original.paths[i].height);
        REQUIRE(loop.paths[i].mm3_per_mm == original.paths[i].mm3_per_mm);
    }
    check_closed(loop);
}

TEST_CASE("AeroWing joins hole contours once and chains holes aligned along Y", "[AeroWing]")
{
    ExtrusionLoop outer = make_loop(rectangle(60, 40));
    std::vector<ExtrusionLoop> holes;
    auto add_hole = [&](double x, double y) {
        Polygon polygon = rectangle(4, 4);
        polygon.translate(Point(scale_(x), scale_(y)));
        polygon.make_clockwise();
        holes.push_back(make_loop(polygon));
        holes.back().set_loop_role(elrPerimeterHole);
    };
    add_hole(10, 8);
    double total_gap = 40. - 12.;
    SECTION("one hole") {}
    SECTION("separate X positions") { add_hole(40, 15); total_gap += 40. - 19.; }
    SECTION("three holes in the same Y direction") {
        add_hole(10, 20);
        add_hole(10, 30);
        total_gap = 6. + 6. + 8.;
        std::reverse(holes.begin(), holes.end()); // independent of input order
    }
    SECTION("reversed exterior") { outer.reverse(); }
    SECTION("reversed hole input") { holes.front().reverse(); }
    SECTION("dense vertices and path boundaries at both attachments") {
        for (ExtrusionLoop *loop : {&outer, &holes.front()}) {
            const ExtrusionPath prototype = loop->paths.front();
            const Points points = prototype.polyline.points;
            const coord_t top = get_extents(loop->polygon()).max.y();
            loop->paths.clear();
            for (size_t i = 1; i < points.size(); ++i) {
                Points edge {points[i - 1]};
                if (points[i - 1].y() == top && points[i].y() == top) {
                    for (double x : {12.1, 12.2, 12.3}) edge.emplace_back(coord_t(scale_(x)), top);
                    if (points[i - 1].x() > points[i].x()) std::reverse(edge.begin() + 1, edge.end());
                }
                edge.push_back(points[i]);
                for (size_t j = 1; j < edge.size(); ++j) {
                    ExtrusionPath path = prototype;
                    path.polyline = Polyline {edge[j - 1], edge[j]};
                    loop->paths.push_back(std::move(path));
                }
            }
        }
    }
    SECTION("segmented hole with different flow") {
        const ExtrusionPath prototype = holes.front().paths.front();
        const Points points = prototype.polyline.points;
        holes.front().paths.clear();
        for (size_t i = 1; i < points.size(); ++i) {
            ExtrusionPath path = prototype;
            path.mm3_per_mm *= 1. + 0.1 * i;
            path.polyline = Polyline {points[i - 1], points[i]};
            holes.front().paths.push_back(std::move(path));
        }
    }
    ExPolygon material(outer.polygon());
    material.contour.make_counter_clockwise();
    double boundary_length = outer.length();
    for (const ExtrusionLoop &hole : holes) {
        material.holes.push_back(hole.polygon());
        material.holes.back().make_clockwise();
        boundary_length += hole.length();
    }
    const double width = outer.paths.front().width;
    const double spacing = Flow::rounded_rectangle_extrusion_spacing(width, outer.paths.front().height);
    // A large grid pitch isolates the hole connectors from ordinary fins.
    REQUIRE(AeroWing::connect_holes(outer, holes, 1000.));
    REQUIRE(outer.loop_role() & elrAeroWing);
    check_closed(outer);
    size_t downward_connectors = 0;
    Points connector_ports;
    for (const ExtrusionPath &path : outer.paths) {
        const Points &points = path.polyline.points;
        if (points.size() == 4 && points.front().x() == points.back().x() &&
            points[1].x() == points[2].x() && points[1].x() == points.front().x()) {
            REQUIRE(points.front().y() > points.back().y());
            REQUIRE(points.front().y() > points[1].y());
            REQUIRE(points[2].y() > points.back().y());
            connector_ports.push_back(points.front());
            connector_ports.push_back(points.back());
            ++downward_connectors;
        }
    }
    REQUIRE(downward_connectors == holes.size());
    // Apart from closing the layer loop, no point may be revisited: the two
    // ports replace the shared junctions, including between chained holes.
    Points traversal;
    for (const ExtrusionPath &path : outer.paths)
        for (const Point &point : path.polyline.points)
            if (traversal.empty() || traversal.back() != point) traversal.push_back(point);
    REQUIRE(traversal.front() == traversal.back());
    traversal.pop_back();
    std::map<std::pair<coord_t, coord_t>, std::vector<size_t>> visits;
    for (size_t i = 0; i < traversal.size(); ++i)
        visits[{traversal[i].x(), traversal[i].y()}].push_back(i);
    for (const auto &visit : visits) {
        REQUIRE(visit.second.size() == 1);
    }
    for (size_t i = 0; i < traversal.size(); ++i) {
        const Line a(traversal[i], traversal[(i + 1) % traversal.size()]);
        for (size_t j = i + 1; j < traversal.size(); ++j) {
            if (j == i + 1 || (j + 1) % traversal.size() == i) continue;
            const Line b(traversal[j], traversal[(j + 1) % traversal.size()]);
            Point intersection;
            REQUIRE_FALSE(a.intersection(b, &intersection));
        }
    }
    for (const ExtrusionPath &path : outer.paths)
        REQUIRE(diff_pl(path.polyline, offset(material, 2.f)).empty());
    const double connectors_length = 2. * total_gap - holes.size() * 2. * spacing;
    REQUIRE(unscale<double>(outer.length() - boundary_length) == Approx(connectors_length).margin(1e-4));
    Polylines generated;
    outer.collect_polylines(generated);
    for (const ExtrusionLoop &hole : holes) {
        size_t cuts = 0;
        for (const Point &port : connector_ports) {
            bool on_hole = false;
            for (const ExtrusionPath &path : hole.paths) {
                const Points &points = path.polyline.points;
                for (size_t i = 1; i < points.size(); ++i)
                    if (Line(points[i - 1], points[i]).distance_to(port) <= 2.) on_hole = true;
            }
            if (on_hole) ++cuts;
        }
        REQUIRE(cuts >= 1);
        // Each hole is followed once except for the short gaps at its ports.
        Polylines hole_paths;
        hole.collect_polylines(hole_paths);
        double covered_length = 0.;
        for (const Polyline &part : intersection_pl(generated, offset(hole_paths, 2.f)))
            covered_length += part.length();
        REQUIRE(covered_length == Approx(hole.length() - cuts * scale_(spacing)).margin(100.));
    }
}

TEST_CASE("AeroWing leaves the exterior intact when a hole cannot be connected", "[AeroWing]")
{
    ExtrusionLoop outer = make_loop(rectangle(30, 20));
    const ExtrusionLoop original = outer;
    Polygon polygon = rectangle(4, 4);
    polygon.translate(Point(scale_(13.), scale_(15.8)));
    polygon.make_clockwise();
    ExtrusionLoop hole = make_loop(polygon);
    hole.set_loop_role(elrPerimeterHole);
    REQUIRE_FALSE(AeroWing::connect_holes(outer, {hole}, 20.));
    REQUIRE(outer.loop_role() == original.loop_role());
    REQUIRE(outer.paths.size() == original.paths.size());
    REQUIRE(outer.paths.front().polyline.points == original.paths.front().polyline.points);
}

TEST_CASE("AeroWing shortens grid fins at the first hole instead of dropping them", "[AeroWing]")
{
    ExtrusionLoop outer = make_loop(rectangle(60, 40));
    Polygon hole_polygon = rectangle(20, 10);
    hole_polygon.translate(Point(scale_(15.), scale_(10.)));
    const double width = outer.paths.front().width;
    const double lane_spacing = Flow::rounded_rectangle_extrusion_spacing(width, outer.paths.front().height);
    double expected_tip_y = 10. - width;
    bool second_hole = false;
    std::vector<AeroWing::StiffenerGrid> grids {{10., 0.}};
    SECTION("rectangular hole") {}
    SECTION("two stiffener grids with a hole") { grids.push_back({12., 45.}); }
    SECTION("sloping lower hole boundary") {
        hole_polygon = Polygon {{scale_(15.), scale_(10.)}, {scale_(35.), scale_(14.)},
                                {scale_(35.), scale_(24.)}, {scale_(15.), scale_(20.)}};
        expected_tip_y = 10. + (20. - lane_spacing - 15.) * 0.2 - width;
    }
    SECTION("hole intersects only the offset pass") {
        hole_polygon = rectangle(4.8, 10);
        hole_polygon.translate(Point(scale_(15.), scale_(10.)));
    }
    SECTION("nearest of two holes along Y") { second_hole = true; }
    hole_polygon.make_clockwise();
    std::vector<ExtrusionLoop> holes {make_loop(hole_polygon)};
    holes.back().set_loop_role(elrPerimeterHole);
    if (second_hole) {
        Polygon polygon = rectangle(20, 6);
        polygon.translate(Point(scale_(15.), scale_(26.)));
        polygon.make_clockwise();
        holes.push_back(make_loop(polygon));
        holes.back().set_loop_role(elrPerimeterHole);
    }
    ExPolygon material(outer.polygon());
    for (const ExtrusionLoop &hole : holes) material.holes.push_back(hole.polygon());
    REQUIRE(AeroWing::connect_holes(outer, holes, grids));
    check_closed(outer);
    bool found_shortened = false, found_unaffected = false;
    for (const ExtrusionPath &path : outer.paths) {
        REQUIRE(diff_pl(path.polyline, offset(material, 2.f)).empty());
        // Rectangular/polygonal boundaries and connectors have fewer vertices.
        if (path.polyline.points.size() < 17) continue;
        coord_t tip_y = path.first_point().y();
        for (const Point &point : path.polyline.points) tip_y = std::max(tip_y, point.y());
        if (path.last_point().x() == scale_(20.)) {
            found_shortened = true;
            REQUIRE(unscale<double>(tip_y) == Approx(expected_tip_y).margin(2e-6));
        }
        if (path.last_point().x() == scale_(10.)) {
            found_unaffected = true;
            REQUIRE(unscale<double>(tip_y) == Approx(40. - width).margin(2e-6));
        }
    }
    REQUIRE(found_shortened);
    REQUIRE(found_unaffected);
}

TEST_CASE("AeroWing preserves holes during slicing and rejects separate islands", "[AeroWing]")
{
    Model model;
    ModelObject *object = model.add_object();
    TriangleMesh solid = make_cube(30, 20, 3);
    bool islands = false;
    SECTION("negative volume") {
        object->add_volume(solid);
        TriangleMesh hole = make_cube(4, 4, 5);
        hole.translate(13, 8, -1);
        object->add_volume(hole)->set_type(ModelVolumeType::NEGATIVE_VOLUME);
    }
    SECTION("enclosed cavity in the mesh") {
        TriangleMesh hole = make_cube(4, 4, 1);
        hole.translate(13, 8, 1);
        hole.flip_triangles();
        solid.merge(hole);
    }
    SECTION("two islands in the mesh") {
        islands = true;
        TriangleMesh island = make_cube(10, 10, 3);
        island.translate(40, 0, 0);
        solid.merge(island);
    }
    if (object->volumes.empty())
        object->add_volume(solid);
    object->add_instance();
    object->ensure_on_bed();
    DynamicPrintConfig config = test_config();
    config.set_key_value("aero_wing_mode", new ConfigOptionBool(true));
    config.set_key_value("spiral_mode", new ConfigOptionBool(true));
    config.set_key_value("wall_loops", new ConfigOptionInt(1));
    config.set_key_value("bottom_shell_layers", new ConfigOptionInt(1));
    config.set_key_value("bottom_shell_thickness", new ConfigOptionFloat(0.));
    config.set_key_value("detect_overhang_wall", new ConfigOptionBool(false));
    Print print;
    print.set_status_silent();
    print.auto_assign_extruders(object);
    print.apply(model, config);
    if (islands) {
        REQUIRE_THROWS_AS(print.process(), SlicingError);
    } else {
        print.process();
        size_t layers_with_holes = 0;
        for (const Layer *layer : print.objects().front()->layers()) {
            if (layer->id() < 1) continue;
            const LayerRegion *region = layer->regions().front();
            if (region->slices.surfaces.front().expolygon.holes.empty()) continue;
            ++layers_with_holes;
            REQUIRE(region->perimeters.entities.size() == 1);
            const auto *wrapper = dynamic_cast<const ExtrusionEntityCollection *>(region->perimeters.entities.front());
            REQUIRE(wrapper != nullptr);
            // GCode::Region::append must forward the loop, not its collection.
            REQUIRE(wrapper->can_sort());
            REQUIRE(wrapper->entities.size() == 1);
            REQUIRE(dynamic_cast<const ExtrusionLoop *>(wrapper->entities.front()) != nullptr);
            const auto perimeters = region->perimeters.flatten();
            REQUIRE(perimeters.items_count() == 1);
            const auto *loop = dynamic_cast<const ExtrusionLoop *>(perimeters.entities.front());
            REQUIRE(loop != nullptr);
            REQUIRE(loop->loop_role() & elrAeroWing);
            check_closed(*loop);
        }
        REQUIRE(layers_with_holes > 0);
    }
}

TEST_CASE("AeroWing paths remain a continuous rising extrusion through SpiralVase", "[AeroWing]")
{
    ExtrusionLoop loop = make_loop(rectangle(30, 20));
    REQUIRE(AeroWing::add_stiffener(loop));
    PrintConfig config;
    config.aero_wing_mode.value = true;
    config.spiral_mode.value = true;
    config.spiral_mode_smooth.value = true; // AeroWing must bypass nearest-segment smoothing.
    config.use_relative_e_distances.value = true;
    SpiralVase spiral(config);
    const Point start = loop.first_point();
    std::ostringstream initial;
    initial << "G1 Z0.2\nG1 X" << unscale<double>(start.x()) << " Y" << unscale<double>(start.y()) << "\n";
    spiral.process_layer(initial.str(), false);
    spiral.enable(true);
    GCodeReader reader;
    reader.apply_config(config);
    reader.parse_buffer(initial.str());
    for (int layer = 2; layer <= 4; ++layer) {
        std::ostringstream gcode;
        gcode << std::fixed << std::setprecision(6) << "G1 Z" << 0.2 * layer << '\n';
        Points targets;
        for (const ExtrusionPath &path : loop.paths)
            for (size_t i = 1; i < path.polyline.points.size(); ++i) {
                const Point &p = path.polyline.points[i];
                targets.push_back(p);
                gcode << "G1 X" << unscale<double>(p.x()) << " Y" << unscale<double>(p.y()) << " E0.1\n";
            }
        size_t moves = 0;
        const std::string output = spiral.process_layer(gcode.str(), false);
        reader.parse_buffer(output, [&](GCodeReader &r, const GCodeReader::GCodeLine &line) {
            if (!line.cmd_is("G1") || line.dist_XY(r) == 0.) return;
            REQUIRE(line.extruding(r));
            REQUIRE(line.has_z());
            REQUIRE(line.new_Z(r) >= r.z());
            REQUIRE(moves < targets.size());
            REQUIRE(line.new_X(r) == Approx(unscale<double>(targets[moves].x())).margin(0.00001));
            REQUIRE(line.new_Y(r) == Approx(unscale<double>(targets[moves].y())).margin(0.00001));
            ++moves;
        });
        REQUIRE(moves == targets.size());
        REQUIRE(reader.z() == Approx(0.2 * layer));
    }
}

TEST_CASE("AeroWing integrates with sliced cube perimeters and invalidation", "[AeroWing]")
{
    Model model;
    ModelObject *object = model.add_object();
    object->add_volume(make_cube(30, 20, 3));
    size_t hole_count = 0;
    SECTION("solid section") {}
    SECTION("one through hole") { hole_count = 1; }
    SECTION("two through holes aligned along Y") { hole_count = 2; }
    for (size_t i = 0; i < hole_count; ++i) {
        TriangleMesh hole = make_cube(4, 4, 5);
        hole.translate(8, 4 + 8 * i, -1);
        object->add_volume(hole)->set_type(ModelVolumeType::NEGATIVE_VOLUME);
    }
    object->add_instance();
    object->ensure_on_bed();
    DynamicPrintConfig config = test_config();
    config.set_key_value("spiral_mode", new ConfigOptionBool(true));
    config.set_key_value("aero_wing_mode", new ConfigOptionBool(true));
    config.set_key_value("aero_wing_stiffener_spacing", new ConfigOptionFloat(10.));
    config.set_key_value("wall_loops", new ConfigOptionInt(1));
    config.set_key_value("bottom_shell_layers", new ConfigOptionInt(2));
    config.set_key_value("bottom_shell_thickness", new ConfigOptionFloat(0.));
    config.set_key_value("top_shell_layers", new ConfigOptionInt(0));
    config.set_key_value("sparse_infill_density", new ConfigOptionPercent(0.));
    config.set_key_value("detect_overhang_wall", new ConfigOptionBool(false));
    config.set_key_value("use_relative_e_distances", new ConfigOptionBool(true));
    config.set_key_value("enable_prime_tower", new ConfigOptionBool(false));
    config.set_key_value("brim_type", new ConfigOptionEnum<BrimType>(btNoBrim));
    config.set_key_value("skirt_loops", new ConfigOptionInt(0));
    config.set_key_value("skirt_height", new ConfigOptionInt(0));
    Print print;
    print.set_status_silent();
    print.set_no_check_flag(true);
    print.auto_assign_extruders(object);
    print.apply(model, config);
    REQUIRE(print.objects().size() == 1);
    print.process();
    size_t fins = 0;
    for (const Layer *layer : print.objects().front()->layers()) {
        const auto perimeters = layer->regions().front()->perimeters.flatten();
        if (layer->id() < 2) {
            REQUIRE(perimeters.items_count() == hole_count + 1);
            for (const ExtrusionEntity *entity : perimeters.entities) {
                const auto *loop = dynamic_cast<const ExtrusionLoop *>(entity);
                REQUIRE(loop != nullptr);
                REQUIRE_FALSE(loop->loop_role() & elrAeroWing);
                check_closed(*loop);
            }
            continue;
        }
        REQUIRE(perimeters.items_count() == 1);
        const auto *loop = dynamic_cast<const ExtrusionLoop *>(perimeters.entities.front());
        REQUIRE(loop != nullptr);
        check_closed(*loop);
        if (layer->id() >= 2) {
            REQUIRE(loop->loop_role() & elrAeroWing);
            REQUIRE(std::count_if(loop->paths.begin(), loop->paths.end(), [](const ExtrusionPath &path) {
                return path.polyline.points.size() >= 17;
            }) == 3);
            ++fins;
        }
        else REQUIRE_FALSE(loop->loop_role() & elrAeroWing);
    }
    REQUIRE(fins > 2);

    // Exercise the real exporter, including loop orientation, seam placement,
    // simplification, and the serial SpiralVase pipeline.
    const auto filename = boost::filesystem::temp_directory_path() / boost::filesystem::unique_path("aerowing-%%%%-%%%%.gcode");
    GCodeProcessorResult result;
    print.export_gcode(filename.string(), &result, nullptr);
    std::ifstream file(filename.string());
    REQUIRE(file.good());
    const std::string exported((std::istreambuf_iterator<char>(file)), std::istreambuf_iterator<char>());
    file.close();
    boost::filesystem::remove(filename);
    GCodeReader reader;
    reader.apply_config(config);
    const std::string layer_tag = ";" + GCodeProcessor::reserved_tag(GCodeProcessor::ETags::Layer_Change);
    int layer_id = -1;
    size_t spiral_moves = 0;
    reader.parse_buffer(exported, [&](GCodeReader &r, const GCodeReader::GCodeLine &line) {
        if (line.raw() == layer_tag) ++layer_id;
        // Skip the bottom, transition, and last layer (which contains ramp-out).
        if (layer_id <= 2 || layer_id + 1 >= int(print.objects().front()->layers().size())) return;
        if (line.cmd_is("G0") || line.cmd_is("G1")) {
            REQUIRE_FALSE(line.retracting(r));
            if (line.dist_XY(r) > 0.) {
                REQUIRE(line.extruding(r));
                REQUIRE(line.has_z());
                REQUIRE(line.new_Z(r) >= r.z());
                ++spiral_moves;
            }
        }
    });
    REQUIRE(spiral_moves > 20);

    // An angle change must regenerate the walls, anchored after the solid bottom.
    const double tilt = 10.;
    config.set_key_value("aero_wing_stiffener_angle", new ConfigOptionFloat(tilt));
    print.apply(model, config);
    REQUIRE_FALSE(print.objects().front()->is_step_done(posPerimeters));
    print.process();
    double first_spiral_z = 0.;
    coord_t first_root_x = 0;
    for (const Layer *layer : print.objects().front()->layers()) {
        if (layer->id() < 2) continue;
        const auto perimeters = layer->regions().front()->perimeters.flatten();
        const auto *loop = dynamic_cast<const ExtrusionLoop *>(perimeters.entities.front());
        REQUIRE(loop != nullptr);
        REQUIRE(loop->loop_role() & elrAeroWing);
        check_closed(*loop);
        const Point &root = last_fin(*loop).last_point();
        if (layer->id() == 2) {
            first_spiral_z = layer->slice_z;
            first_root_x = root.x();
        }
        REQUIRE(unscale<double>(root.x() - first_root_x) ==
                Approx((layer->slice_z - first_spiral_z) * std::tan(tilt * PI / 180.)).margin(2e-6));
    }

    config.set_key_value("aero_wing_stiffener_spacing", new ConfigOptionFloat(20.));
    print.apply(model, config);
    REQUIRE_FALSE(print.objects().front()->is_step_done(posPerimeters));
    print.process();
    for (const Layer *layer : print.objects().front()->layers()) {
        if (layer->id() < 2) continue;
        const auto perimeters = layer->regions().front()->perimeters.flatten();
        const auto *loop = dynamic_cast<const ExtrusionLoop *>(perimeters.entities.front());
        REQUIRE(loop != nullptr);
        REQUIRE(std::count_if(loop->paths.begin(), loop->paths.end(), [](const ExtrusionPath &path) {
            return path.polyline.points.size() >= 17;
        }) == 1);
    }

    // Every second-set setting must invalidate and regenerate the perimeters.
    for (const std::string &key : {"aero_wing_secondary_stiffeners",
                                   "aero_wing_secondary_stiffener_spacing",
                                   "aero_wing_secondary_stiffener_angle"}) {
        if (key == "aero_wing_secondary_stiffeners")
            config.set_key_value(key, new ConfigOptionBool(true));
        else
            config.set_key_value(key, new ConfigOptionFloat(key == "aero_wing_secondary_stiffener_spacing" ? 6. : 30.));
        print.apply(model, config);
        REQUIRE_FALSE(print.objects().front()->is_step_done(posPerimeters));
        print.process();
    }
    bool added_fins = false;
    for (const Layer *layer : print.objects().front()->layers()) {
        if (layer->id() < 2) continue;
        const auto perimeters = layer->regions().front()->perimeters.flatten();
        REQUIRE(perimeters.items_count() == 1);
        const auto *loop = dynamic_cast<const ExtrusionLoop *>(perimeters.entities.front());
        REQUIRE(loop != nullptr);
        check_closed(*loop);
        added_fins |= std::count_if(loop->paths.begin(), loop->paths.end(), [](const ExtrusionPath &path) {
            return path.polyline.points.size() >= 17;
        }) > 1;
    }
    REQUIRE(added_fins);
    config.set_key_value("aero_wing_secondary_stiffeners", new ConfigOptionBool(false));
    print.apply(model, config);
    REQUIRE_FALSE(print.objects().front()->is_step_done(posPerimeters));
    print.process();
    for (const Layer *layer : print.objects().front()->layers()) {
        if (layer->id() < 2) continue;
        const auto perimeters = layer->regions().front()->perimeters.flatten();
        const auto *loop = dynamic_cast<const ExtrusionLoop *>(perimeters.entities.front());
        REQUIRE(loop != nullptr);
        REQUIRE(std::count_if(loop->paths.begin(), loop->paths.end(), [](const ExtrusionPath &path) {
            return path.polyline.points.size() >= 17;
        }) == 1);
    }

    config.set_key_value("aero_wing_mode", new ConfigOptionBool(false));
    print.apply(model, config);
    REQUIRE_FALSE(print.objects().front()->is_step_done(posSlice));
    print.process();
    for (const Layer *layer : print.objects().front()->layers()) {
        const auto perimeters = layer->regions().front()->perimeters.flatten();
        const auto *loop = dynamic_cast<const ExtrusionLoop *>(perimeters.entities.front());
        REQUIRE(loop != nullptr);
        REQUIRE_FALSE(loop->loop_role() & elrAeroWing);
    }
}
