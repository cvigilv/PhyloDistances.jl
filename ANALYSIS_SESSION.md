# Session Handoff — 2026-08-17

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-004: split-extraction

`src/splits.jl` extracts the bipartitions a tree induces on its taxa. A split is a
`BitVector` over `TaxonIndex` positions, oriented so the first taxon is never a member;
`Splits` pairs the masks with the branch length supporting each. A single bottom-up walk
records the branch above every node except the root, canonicalizes each mask, and
accumulates lengths.

`length(symdiff(s1, s2))` is the Robinson-Foulds distance, so the next metric chunk is
mostly assembly.

## Key decisions made

- **Canonical orientation makes the representation inherently unrooted**, which is what
  discharges the promise made by the rooting warning. A rooted tree's two root branches
  describe the *same* bipartition, so they collapse to one split whose length is their
  **sum** — reconstructing the branch an unrooted reading sees. Summing is the general rule
  for any two branches inducing the same bipartition, which also covers degree-2 nodes.
- **A rooted binary tree therefore gives n−3 splits, not n−2.** Rooted metrics
  (Kendall-Colijn) must take their information from path lengths rather than splits.
- **Trivial splits are excluded by default.** They correspond to pendant branches, which
  every tree on the same taxa has, so they distinguish nothing. `trivial = true` keeps them,
  which branch-score-style metrics need since those sum over all branches.
- **Masks are sorted**, so iteration order is a function of the split set alone rather than
  of traversal order. Two trees sharing an index iterate their common splits identically.
  This was found empirically — before sorting, the same split set came out in different
  orders from differently-shaped trees.
- **Set operations throw when the two `Splits` came from different taxon orderings.** Masks
  from different orderings index different taxa, so comparing them would silently answer a
  question nobody asked.
- **`branchlength(node)` is the extension point** for non-NewickTree node types, alongside
  `taxonlabel`. Lengths are `NaN` where the tree records none, which propagates visibly.

## State of the codebase

- Files created: `src/splits.jl`, `test/test_splits.jl`
- Files modified: `src/PhyloDistances.jl`, `test/runtests.jl`
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — 171 tests
- Entry points: no analysis entry point yet.
  `julia --project=. -e 'using Pkg; Pkg.test()'` runs the suite.
- Known issues: none
- Verified during development on realistic primate trees: hand-checked that canonical
  orientation flips the expected splits, that the rooted and unrooted forms of one tree give
  the same split set, and that the n−3 invariant holds up to 200 taxa. Not committed as
  tests; the synthetic fixtures cover the same paths portably.

## Next chunk

CHUNK-005: random-tree-generation

Seeded generators for test fixtures: uniform random topology on n tips, random branch
lengths, and controlled perturbation (apply k random NNI moves, returning the perturbed
tree). The perturbation generator is the main tool for validating metrics with no closed
form — distance should rise with k in expectation.

This belongs in `src/` rather than `test/`, documented and tested, so the benchmark and
comparison scripts can use it too. Verify that generated trees are valid (correct tip
count, every internal node has ≥2 children, no repeated labels) and reproducible under a
fixed seed.

CHUNK-006 (robinson-foulds) is the chunk after. The distance itself is
`length(symdiff(s1, s2))`, already verified against TreeDist; the work is the normalizer
(`n1 + n2`), the weighted variant, and `Distances.result_type`.

Four further metric chunks were added after auditing TreeDist's exports — CHUNK-028
(transfer distance), CHUNK-029 (SPR), CHUNK-030 (hierarchical mutual information),
CHUNK-031 (consensus information). Chunk-ID order no longer matches execution order past
CHUNK-020; follow the `Depends on` fields.

## Reference implementation is now available

TreeDist is the reference for every metric's default convention, and **its source must be
read before implementing one**. Both forms are set up:

- **Source**: `git clone --depth 1 https://github.com/ms609/TreeDist.git` into a scratch
  directory. Formulations in `R/tree_distance_*.R`, heavy lifting in `src/*.cpp`.
- **Installed CRAN build 2.14.1**, callable from `Rscript`. Compiling it against the Nix R
  needed gettext, zlib and libuv paths supplied via a scratch `R_MAKEVARS_USER` Makevars.

Split extraction has been cross-checked against it: seven hand-picked pairs agree with
`TreeDist::RobinsonFoulds`, including a rooted tree against its unrooted twin (RF = 0) and
a polytomy against a resolved tree. Split counts match `TreeTools::as.Splits`.

## Normalization was widened

`normalize` now mirrors TreeDist's `how` rather than being a flag: `false` (raw), `true`
(the metric's own scheme), a **function** combining the two trees' `normalizerinfo`, or a
**real number** used directly as the divisor. A metric supplies
`normalizerinfo(m, convention, tree)` — what one tree carries — and `normalizer` sums the
two by default, which is what TreeDist's `+` combiner does. `normalize = max` / `min` then
express a result relative to the more or less informative tree.

Verified against the reference with a Robinson-Foulds-shaped metric: `normalize = true`
reproduces `TreeDist::RobinsonFoulds(..., normalize = TRUE)` exactly on four pairs
(1, 0.5, 1, 0.6), including the polytomy case where `2(n-3)` would be wrong.

Metrics implement their formula **without** `.FloorNumericalNoise`. Measure the divergence
from the reference during validation and add flooring only if it proves material.

## Watch out for

- **`Bool <: Real` in Julia.** The `normalize` dispatch puts `::Bool` ahead of `::Real` so
  `normalize = true` means "the metric's own scheme" while `normalize = 1` means "divide by
  one". Preserve that ordering if the dispatch is touched.
- **Derive normalizers from TreeDist, never from the literature.** TreeDist normalizes
  Robinson-Foulds by `n1 + n2`, the total splits present in the two trees — not by
  `2(n-3)`. The two agree for binary trees and diverge as soon as either has a polytomy.
  Expect the same trap in other metrics.
- **`.FloorNumericalNoise` zeroes reference values** below
  `sqrt(eps) * max(1, treesIndependentInfo)`. Reproducing a formula is not enough to
  reproduce a value.
- **Length arithmetic across rootings is inexact.** Recovering an unrooted branch by
  addition is a floating-point sum: `0.2 + 0.4 == 0.6000000000000001`. A rooted tree and the
  same tree written unrooted give lengths differing in the last ulp, so weighted RF, branch
  score and any other length-comparing metric must assert with `≈`, never `==`.
- **Do not rely on `Splits` iteration order matching tree structure** — masks are sorted.
  Use the set operations or mask lookup.
- **`splits` excludes trivial splits by default.** A metric that sums over *all* branches
  (branch score) must pass `trivial = true`, or it will silently ignore every pendant
  branch.
- **`NewickTree.degree` is graph degree, not child count** (`degree(leaf) == 2`). Use
  `length(AbstractTrees.children(node))`.
- **Do not reintroduce a bare `using AbstractTrees` or `using NewickTree`** — they both
  export `children`, `getroot`, `isroot` and `print_tree`.
- **The critical path is CHUNK-001 → CHUNK-008**, ending in validated Robinson-Foulds and
  quartet distance, needed for immediate use. Nothing past CHUNK-008 is urgent.
- **Dummy metrics in `test/test_interface.jl` use real trees** built by local `star(n)` and
  `rooted(n)` helpers. Metrics declaring `requiresrooted` must be handed `rooted(n)`.
- **`test/fixtures/` and `scripts/` hold `.gitkeep` placeholders.** Delete each when the
  directory gains real content.
- **Julia 1.12 is the declared floor and there is no LTS back-support.**
