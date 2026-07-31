# Boundary extraction from the surfel grid

## What the problem actually is

It is narrower than it was, and being precise about that is what makes the
solution small.

The floor stage raycasts **down** onto the real part for every cell, so heights
are exact. Ramps come out smooth and stair risers stay crisp. **Vertical
staircasing does not exist in this pipeline.**

What is left is the **XZ outline**. Put a rotated rectangle on the floor and the
walkable cells form a jagged staircase around its edge, because cells are
axis-aligned and the wall is not. That jaggedness is what inflates polygon
counts and produces junk portals.

The whole fix is: turn the jagged cell outline into clean straight lines,
**without ever asking a Part anything.** Every step below reads only the surfel
grid.

**Input:** `FloorData` — surfels carrying position, normal, slope, clearance,
and the `"x:z"` index.

## The pipeline

### 1. Separate into layers by connectivity, not height

Two adjacent cells belong to the same layer when their height difference is
under the step tolerance. Label connected components over that relation.

Do not slice by height bands. A balcony over a courtyard would put both in the
same slice and corrupt everything downstream. Connectivity also handles a spiral
ramp passing over itself for free — it is one component that simply never
becomes adjacent to itself.

Everything from here runs per component.

The tolerance must **exceed the rise across one cell on the steepest walkable
slope**: tan(65°) = 2.14, so the default is 2.2. Anything lower shatters every
steep ramp into one layer per cell.

### 2. Euclidean distance transform

Every walkable cell gets `D`, the distance to the nearest non-walkable cell.
Middle of a room, large. Hugging a wall, 1.

It must be **true Euclidean**, not 4- or 8-neighbour stepping — those give a
diamond or a square kernel, so the offset derived from them comes out wrong on
the diagonals, which is exactly where the rotated walls are.

`D` does triple duty: thickness map, erosion test, and the driver of the
narrow-corridor handling in step 6.

### 3. Trace the contour

Walk the outline of the component's cells to get a closed loop of boundary
cells. Trace from **where surfels end**, never from SVO solid voxels — the
octree is deliberately over-conservative and would inflate the boundary by up to
a leaf.

This captures walls and cliff edges identically. A cell next to a wall and a
cell at a rooftop edge are both just "floor stops here". One mechanism, no
special cases.

*Implementation note:* boundary **edges** are chained rather than boundary cells
walked, because the boundary of a set of cells is a closed loop by construction —
no open end to chase, no tolerance involved. The ordered cells the fit needs
fall out of the edge order.

### 4. Greedy line fit

**Here is where the staircase dies.**

Walk the loop maintaining a best-fit line through the cells accepted so far.
After each new cell, measure the **maximum perpendicular distance** from any
accepted cell centre to that line. While it stays under ~1 stud, keep extending.
When it exceeds, close the segment and start fresh from that cell.

Three details that matter:

- **Maximum, not average.** An average lets a shallow corner hide inside a long
  run.
- **Total least squares** (PCA on the cell centres), not ordinary least squares —
  edges can run near-vertical in XZ and OLS blows up there.
- **Include the wall-hugging cells in the fit.** They define where the wall is.
  Offsetting comes later.

A rotated rectangle: the staircase cells all sit within a stud of one straight
line at the true angle, so they collapse into a single segment. The fit finds
the angle without ever being told it.

**Corners are where the fit fails.** They are a byproduct, not a prerequisite.
This is the piece that killed every previous attempt — corner detection was
being solved separately against real geometry, and it inherited every failure
mode of face interpretation.

### 5. Bias the fit inward

After fitting, translate each line inward until no accepted cell centre lies
outward of it.

Without this a fit can sit outward of the true wall and eat into the clearance
margin, so the safety guarantee becomes probabilistic. With it, it is exact.

### 6. Offset inward, graded

Push each line inward along its normal by `clamp(maxD − margin, 0, r)`, where
`r` is the standard agent radius (~1.5).

Grading matters. A hard "skip the offset if `maxD < r`" creates a discontinuity,
so two adjacent polygons straddling the threshold get boundaries offset by 1.5
and by 0 and no longer share clean edges. Grading is continuous, and narrow
corridors automatically keep their walkable area instead of vanishing.

Two guards:

- **Miter limit.** Acute corners make offset lines intersect arbitrarily far out —
  the classic spike. Past a threshold, bevel instead.
- **Never dilate-then-erode.** Closing would bridge any wall thinner than the
  kernel and weld two rooms together, which is the through-wall portal bug
  manufactured on purpose. Erosion alone can only remove walkable cells, so it
  can never invent connectivity. **That is the safety property the whole thing
  rests on.**

### 7. Corners from intersections

Intersect each adjacent pair of **offset** lines. That is the corner position —
sub-cell accurate, sharper than anything the raster could give.

This is also why offsetting lines beats eroding cells up front. Eroding cells
first bevels every convex corner into a small arc, which the segmenter then
reads as two or three short segments instead of one clean corner. More polygons,
which is the opposite of the goal.

### 8. Clean up

- **Merge the loop seam.** Greedy segmentation is order-dependent and the start
  point on a closed loop is arbitrary, so there is a spurious corner wherever the
  walk began. Test the first and last segments for collinearity and merge.
- **Minimum segment length + near-collinear merge**, for curved geometry that
  would otherwise shatter into dozens of tiny pieces.
- Watch the interaction: very short segments make adjacent offset lines nearly
  parallel, so their intersection is numerically unstable. **Merge before
  intersecting.**

### 9. Severance check

Union-find over surfel adjacency, before and after offsetting. Any component
that gets cut off is the signal that the offset just severed the map. Mandatory,
not optional — it is the only detector for that failure.

Report it; never silently repair it. A severed layer is a tuning failure and has
to be visible as one.

## Why this survives real maps

The only geometry contact in the entire chain is the raycast that built the
surfels. Raycasts do not care whether they hit a Block, a Union, a MeshPart, or
four parts interpenetrating. Everything downstream is arithmetic on a grid.

The old approach needed a readable planar face with a usable CFrame — which
unions and meshes do not have, and which interpenetrating parts corrupt. That is
why it gave us hell, and why this does not.

## What it costs

Effectively nothing. Zero raycasts, zero physics queries — a handful of
arithmetic ops per cell, against the two raycasts per cell the floor stage
already pays. The expensive stage is unchanged.

## Development workflow

Serialize the surfel field once, then iterate boundary code against the cached
snapshot. Seconds per run, not 25 minutes.

## Known limits, stated rather than hidden

- A layer that **spirals over itself** puts two heights in one cell key, and the
  2D raster merges them. Connectivity keeps the layers apart correctly, but the
  raster for that one layer is wrong. Does not arise on the current test scene;
  will arise on a helix ramp.
- The offset is measured against a **standard agent radius**. Ground narrower
  than that is preserved (graded offset, not a width floor), but the boundary it
  produces is an agent-radius-informed one, and a very different body may want a
  different bake.
