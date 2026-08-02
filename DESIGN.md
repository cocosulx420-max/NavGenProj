# Boundary extraction from part-aligned local grids

## What the problem actually is

Staircasing is a **disagreement between the sampling grid and the geometry**. Fill
in squares to record where a rotated wall is, and the record comes out as a
jagged run of squares. The wall is straight; the record is not. Everything
downstream reads the record.

There are exactly two cures, and this pipeline uses both, in this order:

1. **Turn the grid so it agrees with the geometry.** Then nothing is lost and
   there is nothing to recover.
2. **Where the grid cannot agree with everything at once, fit the line back.**

An earlier version of this document described only the second. It assumed one
world-aligned raster over the whole map — `FloorData` and its `"x:z"` index —
so *every* rotated part staircased and the fit was load-bearing everywhere. Every
stage was a fit, every fit had a tolerance, and the tolerances fought each other.

Turning the grid first is what makes the fit cheap. On a part-aligned grid a
part's own rim already lies along whole lattice lines and fits with **zero
residual**. What is left for the fit is the genuinely hard case: some *other*
part's footprint crossing this one's lattice at an angle. Same method, a fraction
of the edges.

**Chain:** `SVO → Floor → LocalGrid → Boundary`.
**Input to this stage:** `LocalGrid`'s per-part grids — **not** `FloorData`.
`Floor` still runs underneath, because `LocalGrid` is built from its surfels and
non-block parts reuse them, but the boundary stage never reads it.

## The pipeline

### 1. Turn the grid — `LocalGrid`

Each **Block** part is sampled on its own local axes, using the full frame
including tilt, so a rotated slab is sampled square on its own incline. Its edges
fall on whole cell lines of its own grid.

**Unions, MeshParts and wedges have no meaningful surface axes**, so they get a
world-aligned fallback grid built from the global surfels. Those are the only
places where a part's *own* rim can still staircase, and they are where the
remaining roughness on SmallMap lives.

This also dissolves a problem the world-raster version needed a whole stage for.
One shared raster flattens a balcony and the courtyard beneath it onto the same
`"x:z"` key, so the old step 1 had to grow layers under a 2D-injectivity
constraint to stop one silently overwriting the other. Per-part grids never mix
two parts' surfaces in one index, so **that stage is gone rather than ported.**

A cell that would have been floor but was killed is kept as a `DeadCell` *with
the instance that killed it*. The boundary stage does not use that attribution —
see "never ask a Part anything" below — but it is the correct place to record it
and clearance volumes will want it.

### 2. Trace the contour

Walk the outline of a grid's live cell mask to get a closed loop of boundary
cells, in the grid's own lattice.

Chain boundary **edges** rather than walking boundary cells: the boundary of a
set of cells is a closed loop by construction, so there is no open end to chase
and no tolerance involved. The ordered cells the fit needs fall out of the edge
order.

Trace from where **cells** end, never from SVO solid voxels — the octree is
deliberately conservative and would inflate the outline by up to a leaf.

**A cliff inside a fallback grid still ends the floor.** A world-aligned grid can
hold two surfaces at once (a mesh with a ledge), so membership alone is not the
test: stop at a neighbour that is missing **or more than `stepTol` away in
height**. A block grid is a single plane and never needs this.

**Chaining at a junction is not a free choice.** A cliff between two cells that
are both in the mask emits an edge from *each* side, on the same lattice edge in
opposite directions — the upper rim and the lower rim are two different
boundaries that coincide in plan. The outline therefore contains zero-width
slits whose ends are vertices with four edges leaving them. Taking whichever
candidate comes first there splices one rim onto the other and drops degenerate
slivers out of the walk, which then fail segmentation and get drawn as junk.
Use the standard face traversal of an embedded planar graph: **arriving along an
edge, leave on the next edge clockwise from it.** With the interior kept on the
left, that walks each rim whole and turns a slit around at its tip.

### 3. Greedy line fit — where the remaining staircase dies

Walk the loop maintaining a best-fit line through the boundary **cell centres**
accepted so far. After each new cell, measure the **maximum** perpendicular
distance from any accepted centre to that line. While it stays under `fitTol`
(one cell), keep extending. When it exceeds, close the segment and start fresh
from that cell.

- **Maximum, not average.** An average lets a shallow corner hide inside a long
  run.
- **Total least squares** (PCA on the centres), not ordinary least squares — a
  run can be near-vertical in the grid's frame, where OLS blows up.
- **Orient the fitted direction along the run's own travel.** A principal axis
  has no sign, and everything downstream reads the sign as meaning something: a
  flipped segment gets its "outward" normal pointing into the floor, and step 5
  then biases the wrong way.

**Corners are where the fit fails.** They are a byproduct, not a prerequisite.
Corner detection against real geometry is what killed every earlier attempt.

### 4. A reversal ends a run, and the residual cannot see it

Round the end of a strip narrower than `fitTol` and the walk comes back down the
other side. Every returning cell is within tolerance of the line fitted to the
side it just left, so **the fit never fails and the run swallows the whole
ring**. Test the travel direction as well: a corner turns, only the far end of a
thin strip reverses.

This is not a corner case. SmallMap's staircases are built from individual
`35 × 2` cell step parts — 36 of them, 3.6% of every walkable cell on the map —
and without this test all 36 collapsed and were dropped.

The same guard belongs on the step 8 merge, where the residual test alone would
otherwise weld the two sides of a thin strip together.

### 5. Bias the fit inward

Translate each line inward until no accepted cell centre lies outward of it.
Without this a fit can sit outward of the cells it was fitted to and hand back
ground that is not walkable, so the safety guarantee becomes probabilistic. With
it, it is exact.

Together with fitting to cell **centres**, this is why the polygons sit inside
the walkable cells rather than on the part's true face — measured at ~0.96 studs
inside on average. That is erosion, which is the safe direction.

### 6. When the fit cannot resolve a ring

`fitTol` is one cell because anything the mask can express as straight *is*
straight. The consequence is that **a strip two cells wide can never be
segmented**: its cell centres form a rectangle one cell deep, and a one-cell-deep
rectangle is within tolerance of a line. That is the tolerance meaning what it
says, not a bug.

Discarding such a ring is the one operation here that can make real ground
disappear, so it never happens silently. Two fallbacks, in order:

- **If the ring's cells are exactly the border of their own bounding box, emit
  that box** — four vertices, exact. This reads cell indices only, so it is as
  valid for a Union as for a Block.
- **Otherwise keep the raw lattice outline.** Jagged, but present. An unfitted
  *outer* ring stands as it is, since cell centres already lie inside the
  walkable cells. An unfitted *hole* is pushed out by half a cell, because its
  centres sit half a cell into the obstacle and the hole would otherwise come out
  too small — an obstacle slightly too big is safe, one too small is not.

### 7. Corners from intersections

Intersect each adjacent pair of lines. Sub-cell accurate, sharper than the
lattice could give.

**Miter limit.** An acute corner throws the intersection arbitrarily far out —
the classic spike. Past `miterLimit` from the corner it replaces, bevel across
instead. Near-parallel lines have no usable crossing; fall back to the foot of
the anchor rather than inventing one.

### 8. Clean up, all of it before any corner is intersected

- **Merge the loop seam.** Greedy segmentation is order-dependent and the start
  point on a closed loop is arbitrary, so there is a spurious corner wherever the
  walk began. Test the first and last segments for collinearity and merge.
- **Absorb short and near-collinear runs.** This is also what handles a flight
  built from *separate* parts: nine stacked blocks with flush ends make one
  straight edge crossed by eight seams, each seam putting a one-cell jog in the
  mask that closes a run. Without the merge, eight spurious corners along a
  straight line.
- **Signed, never absolute.** An absolute test calls a 180° reversal "collinear",
  and the two long sides of a thin ring are exactly that.
- Short runs are what make a corner intersection unstable, which is why all of
  this runs **before** step 7.

## The rule this module obeys: never ask a Part anything

No CFrames, no face planes, no sizes, no raycasts, no overlap queries. The
boundary stage reads cell centres and nothing else. The only contact with real
geometry in the whole chain is the raycasting `Floor` and `LocalGrid` already
paid for, and a raycast does not care what it hit.

This was tested against the alternative and the alternative was rejected. A
version of this stage took its lines from the side planes of whichever part
killed each cell — using the `DeadCell` attribution — and reconstructed corners
by crossing those planes. On SmallMap it was *more* accurate: polygons landed on
the true face within 0.01 studs instead of ~0.96. It was still the wrong answer:

- A Union or a MeshPart has no readable planar face. Its bounding box is not its
  shape, so snapping a boundary to that box claims floor where there is air.
- **Interpenetrating parts report faces that are not surfaces.** A block half
  buried in another block still describes four tidy planes, some of which
  correspond to nothing.

Accuracy on a clean test map is not worth reintroducing the failure mode the
design exists to escape.

## Not implemented, and required before this output is usable

The current output is **geometry only**. It traces where floor physically ends,
not where an agent can stand.

- **The graded inward offset.** Push each *wall* line inward by
  `clamp(maxD − margin, 0, agentRadius)`, graded rather than switched so two
  adjacent polygons straddling a threshold do not offset by different amounts and
  stop sharing an edge. A **dropoff is not offset at all** — an agent may walk
  the lip of a ledge, and standing off from every ledge as though it were masonry
  removed 12.4% of SmallMap's cells, from exactly the places worth keeping. This
  needs a wall-vs-dropoff classification, and `maxD` (local ground thickness)
  needs a **true Euclidean** distance transform — 4- or 8-neighbour stepping
  gives a diamond or square kernel and comes out wrong on the diagonals, which is
  where the rotated walls are.
- **Severance check.** Union-find over cell adjacency, before and after
  offsetting. Report it; never silently repair it — a severed layer is a tuning
  failure and has to be visible as one. Annihilation (zero surviving pieces) must
  be reported separately, because `pieces > 1` is false for zero exactly as it is
  for one.
- **Inter-part links.** Regions are per part and not merged or connected, so the
  output is a pile of correct surfaces with no way between them.

**Never dilate-then-erode.** Closing would bridge any wall thinner than the
kernel and weld two rooms together — the through-wall portal bug, manufactured on
purpose. Erosion alone can remove walkable ground but can never invent
connectivity. **That is the safety property the whole thing rests on.**

## What it costs

The boundary stage is arithmetic on cell centres: **0.05 s** for SmallMap's 87
parts and 50,674 cells, against ~2.5 s for `LocalGrid` and its raycasts. The
expensive stage is the sampling, and it is unchanged.

## Development workflow

Serialize the local grids once and iterate boundary code against the cached
snapshot. **Measure in Studio; do not reason from reading the code.** Every
defect found in this stage was found by running it and dumping cell grids, and
two confident hypotheses that "explained" the symptoms were wrong.

## Known limits, stated rather than hidden

- **Fallback grids staircase, and cannot be fixed by turning the grid.** On
  SmallMap, of 71 boundary edges under 1.5 studs, the great majority are on the
  two Union hosts.
- **A fallback grid's index is keyed on world `x:z`, so two surfaces of the same
  Union at the same `x:z` overwrite each other** and one is silently lost. It
  does not arise on SmallMap — 0 collisions across both fallback grids — but the
  previous version of this document dismissed a case as "does not arise on the
  current test scene" and it then arose on 17.9% of the largest layer, so it is
  written down rather than assumed away.
- **Polygons sit ~0.96 studs inside the true face** (min 0.70, max 16.3 where a
  mask covers only part of its face). Erosion, therefore safe, but it is a real
  accuracy cost of fitting to cell centres.
- A part whose mask produces no ring at all is dropped. Currently 1 grid and 1
  cell of 50,674 on SmallMap.
