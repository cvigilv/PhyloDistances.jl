# Session Handoff — 2026-08-17

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-006: robinson-foulds

`src/robinsonfoulds.jl` adds `RobinsonFoulds` and `WeightedRobinsonFoulds`, the first real
metrics. The weighted one sums, over every split either tree contains, the difference between
its lengths, treating an absent split as length zero.

The topological distance uses **Day's (1985) cluster-table algorithm** in
`src/clustertable.jl`, the same one TreeDist uses, rather than the symmetric difference of
split sets it started as. It is now **4–25× faster than TreeDist** on a single pair; n=1000
went from 26.7 ms to 120 µs.

**`RobinsonFoulds` is validated against TreeDist 2.14.1**: nine hand-picked pairs embedded in
the tests with provenance, plus randomized cross-checks over 1,500 pairs — half with rooted
inputs — with **zero mismatches** on raw and normalized values, and separately against the
split-set formulation it replaced.

## Key decisions made

- **`normalize = true` divides by `n1 + n2`**, the two trees' combined split count, matching
  the reference. Two star trees give `0/0`; TreeDist returns `NaN` and so do we, pinned by a
  test.
- **TreeDist implements no branch-length-weighted RF.** Its family is `RobinsonFoulds`,
  `InfoRobinsonFoulds` (information-weighted, which belongs with the information metrics)
  and `JaccardRobinsonFoulds`; branch-length weighting is phangorn's `wRF.dist`. So
  `WeightedRobinsonFoulds` has no reference, both conventions coincide, and its tests are
  hand-computed. **The same will apply to the Kuhner-Felsenstein branch score.**
- **Weighted RF sums over the union, not the symmetric difference.** The chunk description
  said otherwise, but summing only over the symmetric difference would ignore length
  disagreement on shared splits, which is not the standard definition. Trivial splits are
  included, so pendant branch lengths count.
- **A tree without branch lengths is rejected** by the weighted metric rather than silently
  compared as `NaN`.
- **`result_type` is `Int` for unnormalized RF and `Float64` otherwise**, so an all-pairs
  matrix of raw RF distances comes back as integers.

## State of the codebase

- Files created: `src/robinsonfoulds.jl`, `src/clustertable.jl`,
  `test/test_robinsonfoulds.jl`, `test/test_clustertable.jl`, `benchmark/`
- Files modified: `src/PhyloDistances.jl`, `src/taxa.jl`, `test/runtests.jl`,
  `test/Project.toml`
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — 900 tests
- Entry points: `RobinsonFoulds()(t1, t2)`, `pairwise(RobinsonFoulds(), trees)`
- Known issues: none

## Usable right now

```julia
using PhyloDistances
t1 = readnw("(((Human,Chimp),Gorilla),Orang,Gibbon);")
t2 = readnw("(((Human,Gorilla),Chimp),Gibbon,Orang);")

RobinsonFoulds()(t1, t2)                     # raw distance
RobinsonFoulds(; normalize = true)(t1, t2)   # on [0, 1]
pairwise(RobinsonFoulds(), [t1, t2, ...])    # all-pairs matrix
```

## Reference implementation

TreeDist is the reference for every metric's default convention, and **its source must be
read before implementing one**. Both forms are set up:

- **Source**: `git clone --depth 1 https://github.com/ms609/TreeDist.git` into a scratch
  directory. Formulations in `R/tree_distance_*.R`, heavy lifting in `src/*.cpp`.
- **Installed CRAN build 2.14.1**, callable from `Rscript`. Compiling it against the Nix R
  needed gettext, zlib and libuv paths supplied via a scratch `R_MAKEVARS_USER` Makevars.

The randomized cross-check pattern is worth repeating for each metric: generate trees in
Julia, write Newick and the Julia value to a TSV, and have an R script recompute and
compare. It is not committed, since it needs R.

## Next chunk

CHUNK-007: quartet-distance-exact

The number of four-taxon subsets whose induced topology differs between the trees, by
direct O(n⁴) enumeration over taxon quadruples. This is deliberately the naive algorithm:
it is the correctness oracle for the faster version later, and fast enough for immediate
use at modest tip counts. The docstring must state the complexity plainly.

Verify: identical trees → 0; a resolved tree against a star gives the total quartet count
`C(n,4)`; hand-computed cases on 5 and 6 tips; symmetry. TreeDist does **not** export a
quartet distance — that is the Quartet package — so check whether a reference is available
before relying on one.

After that, CHUNK-008 (validation of RF and quartet) closes the critical path.

## Values are identical to TreeDist, not merely close

`validation/crosscheck.jl` compares the two implementations on 2,140 deliberately awkward
cases — star trees, caterpillars, polytomies at three severities on one and both sides,
rooted against unrooted, identical against maximally perturbed — with **no tolerance**:
integers exactly, floats bitwise, `NaN` meeting `NaN`. **Zero mismatches.** It exits non-zero
if anything differs, so it can gate a release, and should be re-run for every metric.

**Compare in the direction R → Julia.** R's `as.numeric` does not reliably round-trip its own
`%.17g` output — it lands half an ulp away on values such as `92/94`. An earlier harness fed
Julia's output through it and reported two differences that did not exist; it would equally
have hidden real ones. R writes, Julia parses and compares.

## Performance discipline established

Metrics use the most efficient algorithm available, taking TreeDist's choice as the starting
point. Robinson-Foulds is Day's (1985) cluster-table algorithm, in `src/clustertable.jl`, and
is now **4–25× faster than TreeDist** on a single pair.

Three findings that generalize to every later metric:

- **Iterating `AbstractTrees.Leaves` allocates a cursor per step** — 39,700 per comparison at
  200 taxa. Use one typed, iterative walk that collects structure and labels together.
- **Do not build a canonical sorted ordering the result does not depend on.** Sorting labels
  was 64% of the remaining time at 1000 taxa. `TaxonIndex` is needed only where split masks
  must be reproducible; comparison-only numbering suffices otherwise.
- **`@inbounds` was assessed and rejected on evidence.** Profiling blamed string hashing and
  allocation, not array indexing. Do not add it without a profile that blames bounds
  checking. `@threads` is worthwhile at the all-pairs level, not inside one comparison.

Profile with the `profile-performance` skill. Note the Julia MCP server was unavailable this
session, so measurements ran as one-shot scripts through Bash rather than in a persistent
Revise session.

## Watch out for

- **Rooting a tree at a taxon turns its original root into a single-branch pass-through
  node**, whose cluster duplicates its only child's. Counting both inflated every distance
  involving a rooted input — a real bug the hand-picked reference pairs caught but 2,000
  random *unrooted* pairs did not. Any traversal enumerating clusters must skip nodes with
  fewer than two branches below them, and randomized checks must include rooted trees.
- **Check whether TreeDist implements a metric before assuming a reference exists.** It has
  no weighted RF and, as far as the export list shows, no quartet distance either.
- **Derive normalizers from TreeDist, never from the literature.** It normalizes RF by
  `n1 + n2`, not `2(n-3)`; they agree for binary trees and diverge once either has a
  polytomy.
- **`.FloorNumericalNoise` zeroes reference values** below
  `sqrt(eps) * max(1, treesIndependentInfo)`. Metrics here implement the plain formula;
  measure the divergence during validation and add flooring only if material.
- **`Bool <: Real` in Julia.** The `normalize` dispatch puts `::Bool` ahead of `::Real`, so
  `normalize = true` is "the metric's own scheme" and `normalize = 1` is "divide by one".
- **Length arithmetic across rootings is inexact** (`0.2 + 0.4 == 0.6000000000000001`), so
  branch-length metrics must assert with `≈`, never `==`.
- **`splits` excludes trivial splits by default.** A metric summing over *all* branches must
  pass `trivial = true`.
- **`NewickTree.Node(id, data)` leaves `parent` and `children` undefined**, not empty. Build
  with `push!(parent, child)`.
- **`===` has no curried form**; `findfirst(===(x), v)` throws, unlike `==`.
- **Do not reintroduce a bare `using AbstractTrees` or `using NewickTree`** — they both
  export `children`, `getroot`, `isroot` and `print_tree`.
- **`test/fixtures/` and `scripts/` hold `.gitkeep` placeholders.** Delete each when the
  directory gains real content — `test/fixtures/` is due content in CHUNK-008.
- **Julia 1.12 is the declared floor and there is no LTS back-support.**
