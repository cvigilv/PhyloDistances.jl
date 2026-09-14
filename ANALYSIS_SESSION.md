# Session handoff, 2026-09-13

## Project maturity target

`package`, PhyloDistances.jl

## What was just completed

CHUNK-015: phylogenetic-information-metrics

Added shared phylogenetic information as a similarity and matching split information as a
distance. Both score every pair of non-trivial splits with table-backed phylogenetic
information formulas, then use the existing dense assignment solver to select the best
one-to-one matching.

## Key decisions made

- SPI follows TreeDist's conservative rule: incompatible split pairs score exactly zero.
- MSI scores the more informative of the two candidate splits induced by matching and
  differing taxon memberships. MSID subtracts twice the optimal matched score from the two
  trees' summed splitwise information.
- SPI normalization uses mean tree information. MSID normalization uses summed tree
  information. Both `:treedist` and `:primary` currently select the same formulas.
- Exact shared splits stay in the full assignment problem. For both SPI and MSI, committed
  fixtures show that fixing exact pairs first can lower the optimum. MCI's exact-pair
  reduction is not valid here.
- No benchmark files were added. This chunk's verification called for analytical tests and
  TreeDist agreement, both of which now pass.

## State of the codebase

- Files created: `src/phylogeneticinformation.jl`,
  `test/test_phylogeneticinformation.jl`.
- Files modified: `src/PhyloDistances.jl`, `test/runtests.jl`,
  `validation/crosscheck.jl`, `validation/report.md`, `ANALYSIS_PLAN.md`, and this handoff.
- Package loads cleanly: yes.
- Test suite passes: yes, 9,337/9,337 with Julia 1.12.6.
- Reference validation passes: yes. Raw and normalized SPI/MSID agree with TreeDist 2.14.1
  across all 1,140 deterministic and seeded random cases, with zero mismatches.
- Entry points: `SharedPhylogeneticInfo()(tree1, tree2)` and
  `MatchingSplitInfoDistance()(tree1, tree2)`.
- Known issues: `validation/README.md` still overstates which comparisons are bitwise. The
  executable cross-check and generated report describe the tolerance correctly.

## Next chunk

CHUNK-020: maximum-agreement-subtree

Implement maximum agreement subtree size and its information-content variant. It will use
split extraction and the phylogenetic information functions now available. Validate small
cases against exhaustive leaf-subset search and compare both forms with TreeDist.

## Watch out for

- Read `R/tree_distance_mast.R` and its C++ implementation before choosing the algorithm or
  defining what the information variant returns.
- Exact-pair preprocessing is scorer-specific. Do not transfer MCI's optimization to SPI,
  MSI, or another metric without an exchange proof and a brute-force check.
- Put `/nix/store/jaqvbj23b52yl0qgcrrb4ysbxdlqlbv5-R-4.4.2-wrapper/bin` first on `PATH`
  before running `validation/crosscheck.jl`.
- Test files share one namespace. New helpers need file-specific names.
