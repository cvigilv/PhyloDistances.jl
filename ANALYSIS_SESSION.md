# Session Handoff — 2026-08-17

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-005: random-tree-generation

`src/random.jl` provides `randomtree(rng, n)`, `nni!(rng, tree)` and
`perturb(rng, tree, k)`. Trees are drawn uniformly from unrooted binary topologies by
repeatedly splitting a uniformly chosen branch; `perturb` returns a copy with `k`
nearest-neighbour interchanges applied. This lives in `src/` rather than `test/` so the
benchmark and comparison scripts can use it.

`perturb` is the main tool for checking metrics with no closed form: distance from the
original rises with the move count.

## Key decisions made

- **Uniform topology sampling by branch splitting.** Each of a k-taxon tree's `2k-3`
  branches is an equally likely attachment point for taxon k+1, which makes every unrooted
  binary topology equally likely. Verified: 200 draws on 6 taxa gave 87 distinct topologies
  of the 105 possible, against ~89 expected under uniformity.
- **One interchange gives RF exactly 2**, with zero variance — it relocates a single
  subtree, so precisely one split is replaced. Asserted in the suite; it is a sharp check on
  split-based metrics generally, not only on the generator.
- **`perturb` copies rather than mutating**, and `k` moves means *at most* `k` interchanges
  apart: moves are drawn independently, so a later one can undo an earlier one. Distance
  rises with `k` in expectation and saturates just below `2(n-3)`.
- **`Random` added as a dependency** (stdlib), to both the package and test projects. Seeded
  `Xoshiro` gives reproducibility; no external RNG package was needed, since the package
  targets a single Julia version.

## State of the codebase

- Files created: `src/random.jl`, `test/test_random.jl`
- Files modified: `src/PhyloDistances.jl`, `test/runtests.jl`, `Project.toml`,
  `test/Project.toml`
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — 394 tests
- Entry points: no analysis entry point yet.
  `julia --project=. -e 'using Pkg; Pkg.test()'` runs the suite.
- Known issues: none

## Correction to an earlier note

An earlier handoff claimed `NewickTree.degree` is graph degree rather than child count,
citing `degree(leaf) == 2`. **That was wrong** — it came from testing an internal node and
calling it a leaf. `degree(n) = isleaf(n) ? 0 : length(n.children)`: it *is* child count,
and returns 0 for a leaf. Either `degree` or `length(AbstractTrees.children(node))` works.

## Reference implementation

TreeDist is the reference for every metric's default convention, and **its source must be
read before implementing one**. Both forms are set up:

- **Source**: `git clone --depth 1 https://github.com/ms609/TreeDist.git` into a scratch
  directory. Formulations in `R/tree_distance_*.R`, heavy lifting in `src/*.cpp`.
- **Installed CRAN build 2.14.1**, callable from `Rscript`. Compiling it against the Nix R
  needed gettext, zlib and libuv paths supplied via a scratch `R_MAKEVARS_USER` Makevars.

Split extraction is already cross-checked against it: seven hand-picked pairs agree with
`TreeDist::RobinsonFoulds`, including a rooted tree against its unrooted twin (RF = 0) and
a polytomy against a resolved tree.

## Next chunk

CHUNK-006: robinson-foulds

Robinson-Foulds and weighted RF. The unnormalized distance is `length(symdiff(s1, s2))`,
already verified against TreeDist. The work is:

- `normalizerinfo` = the tree's split count, so `normalize = true` divides by `n1 + n2`
  (confirmed against the reference, including a polytomy case where `2(n-3)` is wrong);
- weighted RF, summing branch-length differences over the symmetric difference — which needs
  `splits(tree, index; trivial = true)` if pendant branches are to count;
- `Distances.result_type`, which is `Int` unnormalized and `Float64` normalized.

## Watch out for

- **`NewickTree.Node(id, data)` leaves `parent` and `children` undefined**, not empty. Build
  with `push!(parent, child)` or `Node(id, data, parent)`; assigning `node.children` on a
  fresh node raises `UndefRefError`.
- **`===` has no curried form.** `findfirst(===(x), v)` throws
  `ArgumentError: ===: too few arguments`; write `findfirst(y -> y === x, v)`. `==` does
  curry, which makes the asymmetry easy to miss.
- **Derive normalizers from TreeDist, never from the literature.** It normalizes RF by
  `n1 + n2`, the total splits present in the two trees — not `2(n-3)`. They agree for binary
  trees and diverge once either has a polytomy.
- **`.FloorNumericalNoise` zeroes reference values** below
  `sqrt(eps) * max(1, treesIndependentInfo)`. Metrics here implement the plain formula;
  measure the divergence during validation and add flooring only if it proves material.
- **`Bool <: Real` in Julia.** The `normalize` dispatch puts `::Bool` ahead of `::Real`, so
  `normalize = true` means "the metric's own scheme" and `normalize = 1` means "divide by
  one". Preserve that ordering.
- **Length arithmetic across rootings is inexact** (`0.2 + 0.4 == 0.6000000000000001`), so
  weighted RF and branch score must assert with `≈`, never `==`.
- **`splits` excludes trivial splits by default.** A metric summing over *all* branches must
  pass `trivial = true`, or it silently ignores every pendant branch.
- **Do not rely on `Splits` iteration order matching tree structure** — masks are sorted.
- **Do not reintroduce a bare `using AbstractTrees` or `using NewickTree`** — they both
  export `children`, `getroot`, `isroot` and `print_tree`.
- **The critical path is CHUNK-001 → CHUNK-008**, ending in validated Robinson-Foulds and
  quartet distance, needed for immediate use.
- **Dummy metrics in `test/test_interface.jl` use real trees** built by local `star(n)` and
  `rooted(n)` helpers. Metrics declaring `requiresrooted` must be handed `rooted(n)`.
- **`test/fixtures/` and `scripts/` hold `.gitkeep` placeholders.** Delete each when the
  directory gains real content.
- **Julia 1.12 is the declared floor and there is no LTS back-support.**
