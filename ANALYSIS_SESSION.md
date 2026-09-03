# Session handoff, 2026-09-04

## Project maturity target

`package`, PhyloDistances.jl

## What was just completed

CHUNK-036: clustering-information-performance

MCI and CID now remove exact split pairs before constructing the mutual-information matrix
and solving the assignment problem. The remaining score loop uses a `log2(0:n)` lookup table
with marginal terms hoisted by row and column. CID reuses the same table for both tree
entropy sums.

## Key decisions made

- Profiling confirmed both suspected costs. At 1000 taxa, score-matrix construction took
  24.01 ms and assignment took 21.71 ms of a 48.55 ms MCI call. Scoring precomputed
  contingency counts took 13.81 ms, while count arithmetic without logarithms took 0.12 ms.
- Exact pairs are safe to fix before assignment. An exchange with one unmatched partner
  follows from `MI(A, X) <= H(A)`; an exchange with two matched partners follows from the
  triangle inequality for variation of information. On the seeded 1000-taxon pair, this
  reduces the assignment from 997×997 to 193×193.
- Integer-log lookup and exact-match removal were benchmarked independently. On pre-extracted
  1000-taxon splits, lookup alone took 34.68 ms, exact removal alone took 1.72 ms, and both
  took 1.29 ms, compared with 46.83 ms originally.
- Tree entropy and exact-pair entropy use the same lookup-table arithmetic. This keeps raw
  self-MCI equal to tree entropy, normalized self-MCI equal to `1.0`, and self-CID equal to
  `0.0`, all exactly.
- Bounds-check samples remained in the score loop after the two main changes. Applying
  `@inbounds` only to the loop whose matrix axes and count ranges prove every access valid
  cut the reduced 1000-taxon matrix by 11% and the full call by about 3%.
- No broader tree-ingestion or assignment-solver changes were attempted. After this work,
  taxon indexing and split extraction account for about half of the 1000-taxon call and most
  of the 200-taxon call. Those paths are shared with other metrics.

## State of the codebase

- Files modified: `src/clusteringinformation.jl`, `test/test_clusteringinformation.jl`,
  `benchmark/run.jl`, `benchmark/README.md`, `benchmark/results.md`, `ANALYSIS_PLAN.md`, and
  this handoff. `validation/report.md` was regenerated and remains unchanged in content.
- Package loads cleanly: yes.
- Test suite passes: yes, 9,224/9,224 with Julia 1.12.6.
- Reference validation passes: yes. Raw and normalized MCI and CID still match TreeDist
  2.14.1 across all 1,140 deterministic and seeded random cases, with zero mismatches.
- Seeded benchmark passes: yes. Every Julia value agrees with its R reference.
- Entry points: `MutualClusteringInfo()(tree1, tree2)` and
  `ClusteringInfoDistance()(tree1, tree2)`. The public API and formulas are unchanged.
- Final MCI timings at 10, 50, 200 and 1000 taxa are 9.8 µs, 60.9 µs, 271.1 µs and
  2.62 ms. Final CID timings are 9.5 µs, 57.4 µs, 272.6 µs and 2.64 ms.
- At 1000 taxa, MCI is now 2.2× faster than TreeDist and CID is 2.3× faster. Both allocate
  2.89 MB instead of 10.42 MB. MCI remains 1.3× slower than TreeDist at 200 taxa.
- Known issues: `validation/README.md` still overstates which comparisons are bitwise. The
  executable cross-check and generated report describe the tolerance correctly.

## Next chunk

CHUNK-015: phylogenetic-information-metrics

Implement shared phylogenetic information and matching split information distance using the
split-matching framework and `SplitInfoTable`. Read TreeDist's source first, preserve its
incompatible-pair treatment, and validate raw and normalized forms against TreeDist.

## Watch out for

- Put `/nix/store/jaqvbj23b52yl0qgcrrb4ysbxdlqlbv5-R-4.4.2-wrapper/bin` first on `PATH`
  before running `validation/crosscheck.jl`. The benchmark intentionally uses the R 4.6.1
  installation normally found on `PATH`.
- Exact-match removal relies on both the mutual-information bound and the triangle inequality
  for variation of information. Do not copy it to another scorer without proving the
  corresponding exchange argument.
- MCI's table-backed contingency formula differs from the standalone `mutualinformation`
  evaluation only by floating-point rearrangement. Tests compare every matrix cell to the
  standalone formula with a `1e-14` tolerance, and the TreeDist validation uses the existing
  generalized-RF tolerance.
- Test files share one namespace. Give new helper functions file-specific names.
