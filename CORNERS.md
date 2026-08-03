# Finding corners by reading the staircase

This document specifies how corners are found. It replaces the corner handling in
DESIGN.md, which inferred them from where a line fit happened to fail. Everything
else in DESIGN.md — the local grids, the node classes, the safety rules — stands.

## Why corners cannot be a by-product of fitting

Every version of this stage so far discovered corners by accident: grow a line
until its residual fails, and call the failure a corner. That ordering is the
source of essentially every corner defect this project has had.

- The fit is greedy, so it closes a run at the first node that breaks tolerance,
  which is routinely one node *before* the corner.
- The leftover nodes in the crook then become a run of their own — a chamfer
  cutting the corner off.
- A length threshold gets added to delete the chamfers, and the same threshold
  deletes real corners: measured, `cornerSpan = 4` cut 82 corners and moved the
  boundary away from 693 nodes, the worst by 3.71 studs.
- Tuning the threshold trades one fault for the other. There is no setting that
  removes chamfers and keeps corners, because they are the same shape.

A corner is a property of the ground. It exists before any line is drawn and does
not depend on where a walk happened to start. So it must be **found first, as a
fact**, and the fit must be told where the corners are rather than asked to
guess.

**Two passes, in this order:**

1. Find every corner from the node grid alone.
2. Build edges between consecutive corners.

## What the grid actually shows

A wall that is not aligned to the grid it is sampled on leaves a **staircase**: a
run of boundary nodes that mostly advance along one axis and occasionally step
across to the other.

Two axes per run, decided by the run itself and never assumed:

- **travel** — the axis the run advances along
- **step** — the other one

A vertical staircase is the same object with the roles swapped. Nothing in this
document may name a fixed axis.

**The pattern is not a property of the angle alone.** Two walls at the same angle
half a stud apart quantise differently, so they produce different step patterns
while representing identical geometry. An irregular staircase is still a
staircase. **Therefore no rule here reads cadence at all.**

Three attempts did, and each failed: counting steps, then requiring a step to be
bracketed by travel, then measuring a stepping *rate* and ending the run when it
changed. Two of the three measured as outright regressions against hand-marked
corners. They share one premise — that the arrangement of steps identifies the
wall — and that premise is false. Every such rule disagrees with itself across a
sub-cell offset, and the tolerance that hides the disagreement also hides real
turns: one knob, trading recall against precision monotonically.

What identifies a run is not how it steps but **where it is going**.

## The rules

Three, and the first two need no tolerance.

### 1. Travel sign is fixed

A horizontal staircase goes left or right. Never both. A vertical one goes up or
down, never both. A run that reverses along its own travel axis is not a messier
version of the same staircase; it is a different one.

This is exact — a reversal is proof, not evidence — so it is checked before any
threshold gets a say in something already decided.

**It is not a formality.** DESIGN.md step 4 records this same fact discovered
from the other side, as a pass-two bug: rounding the end of a strip narrower than
`fitTol` brings the walk back down the other side, and every returning cell sits
within tolerance of the line it just left, so the residual *never fails*. 36 of
SmallMap's step parts collapsed entirely until a travel-direction test was added.
Rule 3 below is blind to every one of those. It is cheaper to catch here.

### 2. Step sign is fixed

Our staircase steps down; a rim that steps up is not it. Also exact.

This was already in this document, written as though it only applied when
choosing between branches at a junction. It is not junction-specific — it holds
on every move of every run, and the junction section now just refers to it.

### 3. A run may not drift from its own average heading

The only tolerance in pass one, and it exists for exactly one case: a rim that
turns **without reversing** either sign, which is the shallow-angle turn.

Keep a running average heading over the moves accepted so far, anchored at the
running centroid of the run's nodes. While each new node sits close to that
heading, the run continues *however irregularly it steps* — irregular
quantisation is the same wall. When a node sits too far off it, the run has met a
different staircase.

Measured against the average, not against a chord from first node to last: a
chord is pinned to its two endpoints and a shallow bend between them barely moves
it, which is the same reason DESIGN.md uses maximum rather than average residual.

### The trigger is not the corner

All three rules answer one question — *is the run over* — and none of them says
where the corner is. **The corner is the last travel move**, walked back to.

This is the whole difference from the fit-failure corners this document rejects.
Those let the residual failure *be* the corner, which is why they landed a node
early and left a chamfer in the crook. Here the trigger only rings the bell; the
position comes from the travel-axis rule and is exact. On a skim the two look
like the same idea.

### Diagonal moves are judged on both axes

A diagonal advances travel and step at once, so it is tested as both and is never
invisible. This is what put `case2` beyond every step-pattern rule: its rim turns
via a diagonal, so it never registered a step to bracket, and a 74-node run with
zero steps has no cadence to read either. It has a perfectly well-defined
heading, which is what rule 3 reads.

### 45 degrees needs no special case

At 45 degrees the rim advances on both axes every move with both signs constant,
so rules 1 and 2 never fire and the heading is dead straight — one run the whole
way down, which is correct. It is also symmetric, so it does not matter which
axis is called travel.

### A junction is a question, not an answer

**Measured on SmallMap: 0 of 8116 lattice vertices have more than one boundary
edge starting at them, so the trace there does not branch at all.** That is not
grounds for deleting this section. SmallMap is a clean, deliberately simple map
built to prove the concept; it is not evidence about what the boundary looks like
on a large, messy, vertical one.

The mechanism is real and has been measured before. In the world-raster pipeline
a cliff *inside* a traced set emitted a boundary edge from each side, producing
zero-width slits whose ends are lattice vertices with four edges leaving them --
1146 such adjacencies on this same map. Per-part block grids are single planes,
so they cannot produce one. **Fallback grids can**: a Union or MeshPart floor
carrying two surfaces at different heights has a cliff inside its own grid, and
that is the slit. SmallMap has two Union parts. A complex map will have many.

So this stays, and it stays untested until there is a map that exercises it.

A **junction** is a node where the rim offers more than one continuation. It does
not end the run by itself. It means: *more than one thing continues here, go find
out which of them is my staircase, if any.*

Explore each branch and reject it if:

- **it breaks rule 1 or rule 2** — a branch that travels backwards, or steps the
  opposite way, cannot be this staircase whatever else it does. Nothing special
  about junctions here; it is the same test applied to a candidate continuation.
- **it is not actually connected** — a branch adjacent on screen but unreachable
  through walkable nodes is a different rim
- **it drifts** — a branch heading somewhere the run is not going, by rule 3

If some branch survives, the run continues through the junction and it was never
a corner. If none does, the run ends at the last travel move.

This is what makes a lone rectangular building work without a special case. Its
90-degree turn presents a branch, and that branch fails the step-direction test.

### Walk both ways, and let confidence be asymmetric

From the starting node, walk in **both** directions. This is not only fairness
about where the walk began — it is how the travel axis and stepping rate are
learned in the first place, since a single direction has nothing to judge itself
against.

The two directions will not read equally well. Whichever side is cleaner
establishes the axis and rate; that estimate then becomes the **expectation** for
the messier side. Tolerance is therefore earned rather than fixed: a side that is
ambiguous alone becomes readable once the pattern it should be approximating is
known. This is what "maybe three, maybe four, maybe whatever the other side
allows" means in practice.

### The start node may already be a corner

If the two walks immediately disagree — incompatible step directions, or each
side establishing a different travel axis — that is not a failure to read. It is
the reading: the walk started on the corner. A first-class outcome, since a
randomly chosen start will land on one regularly.

## What this is measured against

Recall alone is a trap. An earlier detector scored 46 of 50 hand-marked corners
while marking one boundary node in nine, and 202 of its 692 corners sat in
provably straight ground. Both numbers are always reported:

- **recall** — hand-marked corners found, matched node for node, no distance
  tolerance
- **precision** — corners landing where the surrounding nodes are collinear,
  which needs no ground truth to measure

The bar to beat is 457 corners with 35 in straight ground.

## Rules this must not break

From DESIGN.md, unchanged:

- **Never ask a Part anything.** No CFrames, no face planes, no sizes, no
  raycasts, no overlap queries. Every rule above reads cell data only.
- **Erosion only.** Ground may be given up, never invented.

## Risks to measure, not to reason about

Three things this design could plausibly get wrong. None is settled by argument;
each needs a bake on a map that contains the case.

- **Curves may shatter.** A curve bends continuously, so rule 3 must eventually
  fire on it — the question is whether it fires at a sane interval or every few
  nodes. Note the failure mode is now length-dependent rather than cadence-
  dependent: drift accumulates with distance, so a long gentle curve breaks into
  a few segments and a tight one into many, which is at least the right shape.
  The crescent-shaped Union on SmallMap is the nearest test case; a real map will
  have far worse.
- **Single-cell noise may fake a corner.** One missing or offset cell on an
  otherwise straight wall creates an unbracketed step. Systematic dropouts of
  this kind existed until the leaf-clipped down-ray in `Floor` was fixed, so the
  noise floor is lower than it was, but isolated artefacts have not been counted.
  Telling one-cell noise from a shallow angle change by pattern alone is the
  brittlest part of this design.
- **Junctions are unexercised.** See above: the handling is specified but nothing
  on SmallMap triggers it, so it is unverified code until a messy map runs.

## Open questions

- **Two compatible branches at one junction.** If a rim genuinely forks and both
  forks stay connected and step the same way, which continues the run? Currently
  undecided.
- **Where exactly the corner vertex goes.** The corner *node* is the last travel
  move. Whether the drawn vertex sits on that node's centre or somewhere derived
  from the two runs' axes is a question for pass two.
- **Which reference rule 3 should use.** Perpendicular distance from the run's
  average heading is what is implemented. Average *direction* alone was the other
  candidate. They behave differently at shallow angles — direction is insensitive
  to a slow bend, position oversensitive on a long straight run — so this is
  settled by measurement, not argument.
- **How much work rule 3 is doing.** If it accounts for nearly every corner, the
  exact rules are probably wrong. If it accounts for none, the tolerance is dead
  weight. Worth reporting per-rule rather than only in aggregate.
