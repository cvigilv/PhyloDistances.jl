# Session Handoff — 2026-08-20 (split-information primitives)

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

The user asked to prioritize information-corrected Robinson-Foulds (TreeDist's
`InfoRobinsonFoulds`), which was not yet a distinct chunk in the plan. Before implementing
anything, I cloned TreeDist (`R/tree_distance_rf.R`, `src/tree_distances.cpp`,
`inst/include/TreeDist/mutual_clustering_impl.h`) to find out exactly what it needs. The
useful discovery: despite living in TreeDist's file next to the true generalized-RF family
(Nye/JRF/MCI/CID/SPI), `InfoRobinsonFoulds` matches splits by exact identity — the same
relationship classic Robinson-Foulds already computes via set intersection — not by an
optimal assignment. So it needs the split-information primitives (CHUNK-013) but *not*
the assignment/matching framework (CHUNK-010, already complete). I added CHUNK-033
(info-robinson-foulds) to the plan, re-prioritized it and its one dependency (CHUNK-013)
ahead of CHUNK-012/014/015, and — with the user's go-ahead — implemented CHUNK-013 this
session. CHUNK-033 itself is next, not yet started.

CHUNK-013 (`src/information.jl`) provides the information-theoretic primitives every
information metric in the plan (CHUNK-014, 015, 020, 030, 031, and now 033) is built from:
`log2rooted`/`log2unrooted` (log-space double-factorial tree counts, avoiding overflow),
`splitinfo` (a split's phylogenetic information content — Thorley, Wilkinson & Charleston
1998), `clusteringentropy` (Shannon entropy of a split as a two-class partition — Meilă
2007), and `mutualinformation`/`jointentropy` (Vinh, Epps & Bailey 2010) for a pair of
splits. All formulas were read from TreeDist's source rather than derived from the
literature description alone, and cross-checked numerically against TreeDist's exact
integer sequences for small n via the MCP Julia session before being committed to tests.

## Key decisions made

- **CHUNK-033 does not depend on CHUNK-010.** This is the single most important scoping
  finding this session — see the plan's 2026-08-20 Working knowledge entry for the full
  derivation. Get this wrong and CHUNK-033 would be implemented as an unnecessary,
  much slower assignment-based scorer.
- **No precomputed table for `log2rooted`/`log2unrooted`.** TreeDist tables these up to 64
  tips because they sit in a tight per-pair, per-split loop; nothing in this package calls
  them in a loop yet, so a table was deferred rather than built speculatively (recorded as
  a new Open Questions entry — revisit once CHUNK-033 or the benchmark suite shows the
  `O(n²)` cost of summing over all of a tree's splits actually matters, the same
  naive-first pattern CHUNK-006 followed for Robinson-Foulds itself).
- **`mutualinformation` divides by `n`** to put it on the same per-taxon scale as
  `clusteringentropy`, which is what makes `jointentropy = H₁ + H₂ - MI` correct rather
  than off by a factor of `n`. Verified directly: `mutualinformation(mask, mask) ≈
  clusteringentropy(mask)` for every test mask.

## State of the codebase

- Files created: `src/information.jl`, `test/test_information.jl`.
- Files modified: `src/PhyloDistances.jl` (new `include` and six new `public` names:
  `log2rooted`, `log2unrooted`, `splitinfo`, `clusteringentropy`, `mutualinformation`,
  `jointentropy` — all unexported, extension-point style like `normalizerinfo`),
  `test/runtests.jl` (registers the new test file), `ANALYSIS_PLAN.md` (CHUNK-013 marked
  complete with implementation notes, new CHUNK-033 added, priority-reorder comment above
  CHUNK-009's expansion block, one new Working knowledge entry, one new Open Questions
  entry, session ledger line).
- Package loads cleanly: yes, Julia 1.12.6 via the MCP session.
- Test suite passes: yes — `Pkg.test()` reports 7505/7505 (up from the prior session's
  5451, plus this session's ~2054 new information-primitive tests). The ~400 uncaptured
  rooting warnings from `test_robinsonfoulds.jl` are pre-existing and already tracked in
  Open Questions — unrelated to this session's changes.
- Not yet committed: everything above is staged but not committed, per the project's
  standing preference to show a draft commit message for review first.
- Entry point: none new — these are internal primitives with no user-facing metric yet.
  `using PhyloDistances: splitinfo, clusteringentropy, mutualinformation, jointentropy,
  log2rooted, log2unrooted` reaches them directly if you want to poke at them.
- Known issues: none.

## Next chunk

**CHUNK-033 (info-robinson-foulds)**, depending only on CHUNK-013 (now complete). Per its
plan entry: distance = `SplitwiseInfo(tree1) + SplitwiseInfo(tree2) - 2 * (info content of
the splits shared identically between the two trees)`, where a tree's `SplitwiseInfo` is
the sum of `splitinfo` over its non-trivial splits. Implement as a direct extension of
CHUNK-006's split-intersection logic (`Splits`'s `intersect`, from `src/splits.jl`),
weighted by `splitinfo` from this session — not as a CHUNK-010 `splitmatching` scorer, and
not via CHUNK-006's fast cluster-table path (which computes shared-split *counts*, not
weighted sums — a new pairwise walk over `intersect(splits(t1, idx), splits(t2, idx))` is
needed). Expect a `SplitwiseInfo`-equivalent helper (sum of `splitinfo` over one tree's
non-trivial splits) to be useful both as CHUNK-033's normalizer and as the aggregate term
in the distance itself — TreeDist's own `SplitwiseInfo` serves exactly both roles.

Validate against `TreeDist::InfoRobinsonFoulds` the same way CHUNK-006 validated classic
RF: RCall-based crosscheck in `validation/`, Nix R 4.4.2 on `PATH` (see prior handoffs for
the exact recipe).

## Watch out for

- Everything already flagged in prior handoffs still applies (Nix R path, `readnw` and
  `SubString`, the ~400 uncaptured rooting warnings, `Bool <: Real`, `NewickTree.Node`
  construction, no bare `using AbstractTrees`/`using NewickTree`, the test-file
  shared-namespace risk — give test helpers file-specific names — and the two
  hand-verification-repeats-the-bug lessons from CHUNK-011).
- **TreeDist's `InfoRobinsonFoulds` R wrapper has fast paths (`.FastDistPath`,
  `.FastManyManyPath`) that batch the all-pairs case via `cpp_rf_info_all_pairs` /
  `cpp_splitwise_info_batch`.** Read past them to `InfoRobinsonFouldsSplits` /
  `robinson_foulds_info` (`src/tree_distances.cpp`) for the single-pair definition — the
  fast paths are TreeDist's own CHUNK-024-equivalent optimization, not a different
  formula, but it is easy to get confused reading the R file top-to-bottom and think the
  batched code path is the canonical definition.
- **`.FloorNumericalNoise` still applies to the TreeDist reference values** (per the plan's
  Reference implementation section) — don't expect bitwise agreement on a value near zero
  without accounting for it, the same caveat that applied to every earlier TreeDist
  cross-check.
- The scratch clone of TreeDist used for this session's source reading
  (`/private/tmp/claude-501/.../scratchpad/TreeDist`) is session-scoped and will not
  persist; re-clone with `git clone --depth 1 https://github.com/ms609/TreeDist.git` into
  a fresh scratch directory next session per the plan's Reference implementation section.
