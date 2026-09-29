#ifndef slic3r_AeroWingGenerator_hpp_
#define slic3r_AeroWingGenerator_hpp_

#include "../ExtrusionEntity.hpp"

namespace Slic3r::AeroWing {

// All entry points accept a common XY orientation around Z, in [-180, 180]
// degrees: 0 points along +Y, +90 along -X, -90 along +X. Spacing, tilt,
// grid_origin_x and hole connections use this rotated frame, anchored at the
// object origin. The default preserves the original orientation.
struct StiffenerGrid {
    double spacing_mm;
    double angle_degrees;
};

// Insert one two-pass fin into an existing exterior loop. The attachment is at
// the lowest intersection with a +Y line offset from the XY bounding box centre.
// The preceding exterior segment is replaced by a direct connection to the fin.
// G-code export aligns the loop start/end with the preceding nozzle position.
// The remaining exterior segments retain their flow. Returns true if a fin was
// added, its shifted axis misses the section, or there is insufficient room.
// In both skipped cases the loop is unchanged. Returns false, without modifying
// the loop, for invalid inputs or unsupported geometry (e.g. holes or bridges).
// Angle is around Y, measured from +Z toward +X. Height is in mm above the
// first spiral layer; the attachment shifts by height * tan(angle) along X.
bool add_stiffener(ExtrusionLoop &loop, double angle_degrees = 0., double height_from_first_spiral_layer = 0.,
                   double orientation_degrees = 0.);

// Fill the section with parallel fins on a fixed X grid. Spacing is measured
// along X in mm; grid_origin_x uses scaled object coordinates and stays fixed
// across layers. Each layer shifts the grid by height * tan(angle).
// Fins without enough clearance are skipped. False leaves the entire loop
// unchanged; true may also leave it unchanged if no fin fits.
bool add_stiffeners(ExtrusionLoop &loop, double spacing_mm, double angle_degrees = 0.,
                    double height_from_first_spiral_layer = 0., coord_t grid_origin_x = 0,
                    double orientation_degrees = 0.);

// Combine grids in perimeter traversal order, sharing clearance checks so
// coincident or neighbouring fins are not extruded twice.
bool add_stiffeners(ExtrusionLoop &loop, const std::vector<StiffenerGrid> &grids,
                   double height_from_first_spiral_layer = 0., coord_t grid_origin_x = 0,
                   double orientation_degrees = 0.);

// Join every hole from the +Y side at its central X, chaining downward along
// -Y through preceding holes when needed, then add the
// regular grid where it fits in the material. Each hole is traversed once, with
// a short opening between distinct entry/exit ports to avoid backtracking.
// Hole contours may contain overhang perimeter segments; their roles and flows
// are preserved, including on connections originating from those segments.
// Holes too small for distinct entry/exit ports are omitted and treated as closed.
// Other failures leave the exterior intact.
bool connect_holes(ExtrusionLoop &outer, const std::vector<ExtrusionLoop> &holes,
                   double spacing_mm, double angle_degrees = 0.,
                   double height_from_first_spiral_layer = 0., coord_t grid_origin_x = 0,
                   double orientation_degrees = 0.);

bool connect_holes(ExtrusionLoop &outer, const std::vector<ExtrusionLoop> &holes,
                   const std::vector<StiffenerGrid> &grids,
                   double height_from_first_spiral_layer = 0., coord_t grid_origin_x = 0,
                   double orientation_degrees = 0.);

} // namespace Slic3r::AeroWing

#endif
