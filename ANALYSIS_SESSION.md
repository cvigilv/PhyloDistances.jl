# Session handoff, 2026-09-03

## Project maturity target

`package`, PhyloDistances.jl

## What was just completed

CHUNK-014: clustering-information-metrics

Added `MutualClusteringInfo`, a `TreeSimilarity`, and `ClusteringInfoDistance`, a
`TreeMetric`. Both maximize mutual information over a one-to-one matching of non-trivial
splits. CID uses the same optimum and computes `H₁ + H₂ - 2 MCI`, where each `H` is the
summed clustering entropy of one tree. The benchmark harness now compares both metrics
against TreeDist on the same seeded Newick pairs at 10, 50, 200 and 1000 taxa.

## Key decisions made

- Read TreeDist's `R/tree_distance_info.R`, `src/tree_distances.cpp`, and
  `inst/include/TreeDist/mutual_clustering_impl.h` before implementation. TreeDist derives
  CID from the matching that maximizes MCI rather than solving a separate minimum-cost
  assignment over variation-of-information scores.
- MCI with `normalize = true` divides by mean tree entropy. CID divides by summed tree
  entropy. These are TreeDist's documented defaults.
- The score matrix reuses packed split words from `generalizedrf.jl`, precomputes each
  marginal split size, and computes pairwise intersections without allocating temporary
  masks.
- `SplitInfoTable` is not relevant here. It accelerates phylogenetic split information;
  clustering entropy and mutual information use split counts and contingency tables.
- Equal and complementary bipartitions return their clustering entropy directly from the
  internal count-based mutual-information helper. This makes self-similarity exact and CID
  exactly zero on identical trees without adding numerical flooring.

## State of the codebase

- Files created: `src/clusteringinformation.jl`, `test/test_clusteringinformation.jl`, and
  `benchmark/clusteringinfo.R`.
- Files modified: `src/PhyloDistances.jl`, `src/information.jl`, `test/runtests.jl`,
  `benchmark/run.jl`, `benchmark/README.md`, `benchmark/results.md`, `.gitignore`,
  `validation/crosscheck.jl`, `validation/report.md`, `ANALYSIS_PLAN.md`, and this handoff.
- Package loads cleanly: yes.
- Test suite passes: yes, 7824/7824 with Julia 1.12.6.
- Reference validation passes: yes. Raw and normalized MCI and CID matched TreeDist 2.14.1
  across 1,140 deterministic and seeded random cases through 400 taxa, with zero mismatches.
- Performance benchmark completed: yes, using TreeDist 2.14.1 under R 4.6.1. Every
  benchmarked value agreed. At 1000 taxa, MCI takes 48.66 ms here versus 5.89 ms in
  TreeDist, 8.3× slower; CID takes 48.31 ms versus 7.64 ms, 6.3× slower.
- Entry points: `MutualClusteringInfo()(tree1, tree2)` and
  `ClusteringInfoDistance()(tree1, tree2)`. Both accept `convention` and `normalize` in
  their constructors and work with `pairwise`.
- Known issues: the clustering metrics have a serious performance gap from 200 taxa onward.
  `validation/README.md` also still overstates which comparisons are bitwise. The executable
  cross-check and generated report describe the tolerance correctly.

## Next chunk

CHUNK-036: clustering-information-performance

Profile MCI/CID, then test integer-logarithm lookup and exact-match removal. TreeDist uses
both techniques. Keep only optimizations supported by the profile, and rerun the full
benchmark, test suite and reference cross-check afterward.

## Watch out for

- Reference validation used the Nix R 4.4.2 installation. Put
  `/nix/store/jaqvbj23b52yl0qgcrrb4ysbxdlqlbv5-R-4.4.2-wrapper/bin` first on `PATH` before
  running `validation/crosscheck.jl`. The benchmark used TreeDist 2.14.1 under the R 4.6.1
  found normally on `PATH`; keep that environment fixed when comparing before and after.
- The TreeDist source clone used this session is under `/tmp/treedist-chunk014.JcOxEX` and
  should be treated as temporary.
- Test files share one namespace. Give any new test helper a specific name rather than
  reusing generic helpers from another file.
- TreeDist quantizes assignment costs and floors small information distances. Keep using
  `_closeenough` for reference validation rather than requiring bitwise agreement.
