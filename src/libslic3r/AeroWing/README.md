# AeroWing: parallel stiffeners

Enable **Others > AeroWing > Enable AeroWing mode** and **Special mode > Spiral
vase**. Start with an upright prism (for example a 30 × 20 × 20 mm box),
one material, one instance and a few solid bottom layers. Rebuild with
`cmake --build build --target BambuStudio -j 4` after changing the code.

**Stiffener spacing** sets the distance along X between successive stiffeners
in millimetres. A fixed grid centred on the object's X origin covers the entire
section on every spiral layer. For example, a 100 mm wide section with 20 mm
spacing has five candidate stiffeners, subject to clearance at each attachment.
At the lowest intersection of each grid line with the original exterior path,
the generator inserts an inward fin along +Y. Each connection omits the return
to its attachment point and all trailing exterior vertices closer than the
normal extrusion spacing. Existing fins are protected from this trimming.
The entry is interpolated on the actual contour segment at that distance from
the attachment. The preserved skin follows that segment up to the entry, so a
small concave detail does not turn the connection into an exterior chord when
the positions of polygon vertices change between layers.
Fins are inserted in perimeter traversal order for either winding direction.
**Stiffener angle** (enabled only in AeroWing mode) ranges from -89° to +89°:
it tilts the longitudinal stiffener around Y, with 0° vertical along Z.
Positive angles shift higher layers toward +X and negative angles toward -X.
The grid offset is `(slice_z - first_spiral_slice_z) * tan(angle)`, measured from
the first layer after the solid-bottom layer count and thickness requirements.
Thus all stiffeners are parallel, and the grid remains fixed relative to the
object even if individual section bounds change. New grid lines can enter the
section at upper layers when other lines leave it; coverage is not limited to
the stiffeners intersecting the first layer.
Within each XY layer, the fin still points along +Y. Its length is the distance
to the next boundary along +Y at the shifted X position minus one
extrusion width. If the shifted axis no longer crosses the section, the layer
keeps its ordinary perimeter without a fin or slicing error. The same applies
when the section is too short or narrow, or the fin or its connection would
leave the contour: that fin is skipped. Fins too close to an existing fin are
also skipped. If none fit, the original perimeter is preserved. Changing the
angle or spacing regenerates the perimeters. Each fin has two adjacent passes
at the normal extrusion spacing, a tapered attachment and a rounded tip. For an
object extruded along Z, successive fins form a longitudinal stiffener.

**Enable a second set of stiffeners** adds another grid, with independent
**Second set stiffener spacing** and **Second set stiffener angle** controls.
These controls are enabled when both AeroWing and the second set are enabled.
The second set is off by default; its initial spacing is 20 mm and its angle is
-45°. Both sets use the same object X origin and first spiral layer as their
reference, with spacing along X and tilt around Y. Their candidate fins are
combined in perimeter traversal order into the same continuous loop. Coincident
fins are printed once; where fins from either set are too close, only the first
one in traversal order is kept. Both sets use the same hole and connector
clearance rules. Changing either set regenerates the perimeters.

The exterior and all fins are paths within **one closed ExtrusionLoop**. Bottom
layers use normal slicing. No STL/3MF mesh is edited. Disabling AeroWing
invalidates slicing and restores standard vase behavior.

**Close top layer** caps the final two layers with solid infill using the top
surface settings. Both layers print without spiral motion or stiffeners, and
holes in the model remain open. The cap layer count is fixed at two; no extra
configuration is needed. Existing solid bottom layers are preserved on short
objects.

## Holes

A connected section with N holes normally produces one exterior perimeter and
N hole perimeters. AeroWing joins them into a single extrusion loop. For each
hole, its bounding box centre X defines a connector along Y: depart from
the +Y side of the exterior toward -Y, follow the hole perimeter, and return
on an adjacent pass. Each hole is followed once, except for a short opening
between the connector's entry and exit; its separate extrusion entity is removed.
Original perimeter flow values are preserved.
The exterior is traversed counter-clockwise and holes clockwise. For an arrival
from +Y, the outgoing lane is on the +X side and the return follows the central
axis. Each attachment uses two distinct points separated by the normal extrusion
spacing. All intervening contour vertices are removed, with the cut interpolated
on the boundary segment, just as for ordinary fins. The connection replaces
that short contour portion on both the parent boundary and the hole, avoiding
backtracking and crossings at the junction. The
lane order is never swapped to fit a narrow passage, and connectors that would
intersect an existing connector are rejected.

When another hole blocks the route toward -Y, the connector starts at the lower
boundary of the preceding hole above it. Successive holes are therefore visited on the same connected
branch, with the return pass leading back to the exterior. All connectors are
checked against the material between the contours, so none cross a hole interior.
Their positions follow the hole centres on each layer independently of the angle
and spacing of the ordinary grid. A grid fin meeting a hole is shortened to the
first hole boundary, keeping its rounded return tip and the usual one-extrusion-
width clearance to the perimeter centreline. The limit accounts for both passes,
including sloping hole boundaries and holes touching only the offset pass.
Grid fins that conflict with a dedicated connector or have insufficient room
even after shortening are omitted.

Every hole must be connected. Insufficient clearance for a hole connector causes
an explicit slicing error, leaving the original geometry intact. It does not
silently discard that hole's perimeter. The current route follows -Y from an upper
boundary; it does not detour around an obstruction in X. Separate disconnected
islands remain unsupported. On solid-bottom layers the ordinary hole perimeters
are printed normally, before the continuous spiral traversal begins.

Placement follows the slicer's object coordinates; there is no airfoil or chord
recognition. Spacing is measured along X, not perpendicular to the inclined
stiffener planes.

Unsupported sections (multiple islands or bridge
paths) produce a slicing error rather than silently omitting the fin. AeroWing
preserves ordinary slicing contours so it can connect holes and detect islands
that standard vase slicing would otherwise discard/fill.

Nearest-point seam relocation and path simplification are bypassed for the
generated loop. Spiral XY smoothing is disabled while AeroWing is active, even
if enabled in a saved profile, because it can confuse the two sides of the fin.
The existing SpiralVase Z ramp and bottom/top extrusion transitions remain in
use. Deposition at the junctions requires physical validation. The prototype is
intended for constant sections; adhesion and support on tapered/twisted wings
have not been validated. The direct connection replaces the final portion of
the exterior skin within one extrusion spacing of each attachment.

## Regression tests

```sh
cmake -S . -B build -DSLIC3R_BUILD_TESTS=ON
cmake --build build --target aero_wing_tests -j 4
ctest --test-dir build -R '^aero_wing_tests$' --output-on-failure
```

Tests cover single and multiple holes, holes chained along Y, preservation of
hole perimeters, multiple fins, grid spacing and coverage at height, the direct
connection and remaining skin, winding direction, segmented loops, rejected
geometry, continuous extrusion through the SpiralVase postprocessor and the full
G-code exporter, sliced cube perimeters and invalidation when the angle, spacing
or mode changes.
