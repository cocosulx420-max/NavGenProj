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
while representing identical geometry. Any rule that reads the raw pattern will
disagree with itself across a sub-cell offset. What survives the offset is the
**rate**: a 20-degree wall steps about once every three travels wherever it sits.
Rules are stated in terms of rate, and raw patterns are only ever evidence.

## The rules

### A step must be bracketed by travel

A legitimate staircase reads travel…travel, step, travel…travel, step. A step
that is not followed by a resumption of travel was never part of the staircase,
and the run ended at the **last travel move** — which is the corner.

This replaces counting steps. "Three steps in a row" was only ever a proxy for
this, and a proxy with an arbitrary constant in it. One unfollowed step is the
same fact as three.

### Lattice phase gets the benefit of the doubt

Because of the sub-cell offset above, a single unfollowed step is *suspicion*,
not proof. On meeting one, keep walking a short distance and ask whether the
staircase **resumes with the same character**: same travel direction, same step
direction, and roughly the same stepping rate.

- resumes the same → it was phase, the run continues through it
- rate or direction changed → the wall genuinely turned

The lookahead is not a fixed constant. See "confidence is asymmetric" below.

### 45 degrees needs no special case

At 45 degrees the pattern alternates strictly — travel, step, travel, step — so
every step is bracketed and the run continues the whole way down, which is
correct: a 45-degree rim is one straight staircase. It is also symmetric, so it
does not matter which axis is called travel; both choices give the same answer.
Near-45 rims produce the occasional doubled step from phase alone, which the rate
check absorbs.

### A junction is a question, not an answer

A **junction** is a node where the rim offers more than one continuation. It does
not end the run by itself. It means: *more than one thing continues here, go find
out which of them is my staircase, if any.*

Explore each branch and reject it if:

- **it steps the wrong way** — our staircase steps down, a branch that steps up
  cannot be the same staircase whatever else it does
- **it is not actually connected** — a branch adjacent on screen but unreachable
  through walkable nodes is a different rim
- **it is pure step** — no travel at all, which is the bracketing rule applied to
  a branch

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

## Open questions

- **Two compatible branches at one junction.** If a rim genuinely forks and both
  forks stay connected and step the same way, which continues the run? Currently
  undecided.
- **Where exactly the corner vertex goes.** The corner *node* is the last travel
  move. Whether the drawn vertex sits on that node's centre or somewhere derived
  from the two runs' axes is a question for pass two.
- **Near-45 rims.** The rate check should absorb phase-induced doubled steps.
  This needs measuring on the real map rather than reasoning about.
