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

### 1. Separate into layers by connectivity **and** 2D injectivity

Two adjacent cells belong to the same layer when their height difference is
under the step tolerance. Label connected components over that relation.

Do not slice by height bands. A balcony over a courtyard would put both in the
same slice and corrupt everything downstream.

**Connectivity alone is not enough, and the reason is not exotic.** This document
used to claim a balcony over a courtyard was safe because the two are separate
components. That holds only while they are *disconnected*. Ramp them together —
the ordinary case, and the case this pipeline exists for — and they are one
component, which then gets flattened into one 2D raster where both floors
compete for the same `"x:z"` key. The loser is silently dropped and every stage
downstream reads a floor plan that is part ground and part balcony.

Measured on SmallMap before the fix: **7038 of 39407 cells in the largest layer
held two or more heights** — 17.9%, all of them more than 5 studs apart, 6243
more than 10 apart, worst pair 36.5 studs. The visible result was one
39407-cell layer whose outer ring was an 8-vertex quad the size of the map
footprint.

So a cell may join a layer only if it is adjacent within the step tolerance
**and** the layer does not already occupy that `x:z` at an incompatible height.
A refused cell is not discarded — it seeds the next layer, which is exactly the
balcony peeling off the courtyard. Every layer is now injective on `x:z` by
construction, which is the precondition steps 2–8 always silently assumed.

Where a **ClipRamp** covers a cell it *is* the walkable surface and the steps
beneath it are not; the risers are dropped. Keeping both makes every step its
own micro-layer.

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

**A cliff inside a layer still ends the floor.** Membership is not the test. A
layer is connected in 3D but its 2D projection need not be: a ramp climbing
alongside the floor it departed from is one layer, passes the injectivity test,
and still had cells 1.7 studs apart horizontally sitting 15 studs apart
vertically. The raster read that drop as walkable ground and ran polygons
across it. Floor stops at a neighbour that is missing **or a cliff away**.

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

### 6. Offset inward, graded — **walls only**

**A dropoff is not offset at all.** The offset exists so an agent does not clip
a wall or snag on a corner. A cliff edge presents neither, and walking the lip
of a ledge is legitimate. Standing off from every ledge as though it were
masonry removed 12.4% of SmallMap's cells, and removed them from exactly the
places worth keeping — balcony rims, platform edges, stair heads — where the
ground is already narrowest.

Boundary edges are classified from the surfel field alone, no Part consulted: a
**wall** is a neighbouring column whose surface stands *above* us, a **dropoff**
is one whose surface lies *below*, and open air is a dropoff too. The ambiguous
case falls to wall, because over-eroding is safe and under-eroding is not.

A run is wall or dropoff and never a blend — a mixed segment has no single
correct push — so the greedy fit breaks on the class change and merging refuses
to cross one. The classification is despeckled first: it is per cell and reads a
neighbouring column, so it picks up grit, and each speck forces a break that
merging can then never undo.

Push each **wall** line inward along its normal by `clamp(maxD − margin, 0, r)`,
where `r` is the standard agent radius (~1.5).

`maxD` is the local **thickness of the ground behind the line**, and it must be
measured by marching inward. Reading `D` at the boundary cells themselves is
worthless: a boundary cell is 4-adjacent to a non-walkable cell by construction,
so its `D` is always exactly 1, which pins the push at `clamp(1 − margin, 0, r)`
forever and makes the grading below a constant.

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

Counting the pieces of what survives cannot see a layer that did not survive at
all: annihilation is zero pieces, and `pieces > 1` is false for zero exactly as
it is for one. That is the worst form of the failure and it was the one case
that passed silently, so it is reported separately.

### 10. Inter-layer links

A layer is a 2.5D raster and a raster cannot express "upstairs". Step 1 is what
makes each raster mean anything, but it also cuts the ramp free of the floor it
climbs from, so the output is a pile of correct surfaces with no way between
them.

The cut is recoverable exactly, because the predicate that made it undoes it:
two cells 4-adjacent in XZ and within `stepTol` are walkable neighbours, so if
they landed in different layers that is a layer boundary drawn through
traversable ground — a ramp foot, a stair head, a balcony meeting its walkway.
A drop larger than `stepTol` is **not** a link, which is what keeps cliffs
cliffs.

Crossings are clustered by adjacency, so two staircases between the same pair of
layers stay two links instead of one averaged point in the middle of neither.

### 11. Vertex heights, along a walkable slope

Steps 2–8 are entirely 2D. Heights come back at the very end, one per ring
vertex, from the nearest cell of that ring's own layer.

Nearest-in-XZ alone is wrong, and visibly so. A layer legitimately spans a whole
building — ground, stairs and roof are one connected, `x:z`-injective component —
so a vertex sitting on a roof edge has ground cells a stud away in XZ and ten
studs below, and picks them. On SmallMap **58 of 102 holes had two consecutive
vertices 2–4 studs apart on the ground and 10 studs apart in height**. Drawn,
that is a fence of vertical bars down every facade, and it is the same class of
bug as the old `cells[1].y` fallback rather than a new one.

The constraint is just the walkable slope. A ring is a loop of walkable
boundary, so over `dxz` studs of ground the height may move by at most `stepTol`
per stud — the same 65° limit step 1 uses. Walk the loop and let each vertex
take the nearest cell its predecessor's height can reach.

Two details:

- **Seed on what the ring mostly agrees on**, not on its single closest vertex.
  A hole ring around a building footprint has one vertex whose nearest cell *is*
  the roof, at distance zero; pinning there drags the whole ground-level loop up
  and leaves the drop on the edge the walk closes over. Seed on the largest
  cluster of unconstrained picks within a step of each other.
- **The constraint may never fail a vertex.** Where no cell of the layer is
  reachable, the unconstrained nearest still stands and is counted in
  `stats.heightUnreachable`. A wrong height is a drawing artefact; no height is
  a hole in the output.

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

- Two heights in one cell key is now handled by the injectivity rule in step 1.
  It was listed here as a spiral-ramp case that "does not arise on the current
  test scene"; it arose, on 17.9% of SmallMap's largest layer.
- **Polygon count is up, not down, since dropoff classification landed**:
  324 → 491 segments on SmallMap, bevels 18 → 140. Breaking a run at every
  wall/dropoff transition fragments the outline, and the bevels are the miter
  limit correctly handling real notches where a wall line pushed 1.5 meets a
  ledge line pushed 0. Honest geometry, but more of it than before.
- **Small regions stay isolated**: 22 regions collapse to 17 islands via links.
  The remainder are crates and small platforms that likely need *jump* links,
  not walk links. Connectivity is never invented to improve that number.
- The distance transform still measures thickness **across** a same-layer cliff,
  so the offset there over-estimates available ground and pushes further than it
  needs to. Over-pushing is erosion, which is safe, so it is left for now.
- A staircase with **no ClipRamp over it** still fragments into micro-layers and
  draws picket-fence rings. On SmallMap the stairs near (81, 105) have none.
- The offset is measured against a **standard agent radius**. Ground narrower
  than that is preserved (graded offset, not a width floor), but the boundary it
  produces is an agent-radius-informed one, and a very different body may want a
  different bake.
