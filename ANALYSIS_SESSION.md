# Session Handoff — 2026-08-20 (info-robinson-foulds, benchmarked)

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-033: `InfoRobinsonFoulds`, the information-corrected Robinson-Foulds distance
(TreeDist's `InfoRobinsonFoulds`).

Before writing any code, I re-cloned TreeDist into scratch (the prior session's clone did
not persist) and read `robinson_foulds_info` in `src/tree_distances.cpp` directly, to
confirm the prior session's scoping conclusion: despite `InfoRobinsonFouldsSplits` calling
`GeneralizedRF(...)` — the same wrapper the true assignment-matched metrics (Nye, JRF, MCI,
CID, SPI) use — the C++ underneath is a plain double loop over both trees' splits that
matches by exact identity (or complement) and sums each match's phylogenetic information;
there is no assignment solve. So the metric is a direct extension of classic
Robinson-Foulds's split-intersection logic, not an instance of CHUNK-010's split-matching
framework, exactly as the prior session inferred from the R source alone.

Implementation lives in `src/robinsonfoulds.jl`, next to `RobinsonFoulds` and
`WeightedRobinsonFoulds` (not in `generalizedrf.jl`, which holds the assignment-based
family). `_splitwiseinfo(s::Splits) = sum(splitinfo, s; init = 0.0)` sums a tree's own
non-trivial splits' information content — TreeDist's `SplitwiseInfo`, serving as both the
metric's own scale and its `normalizerinfo`. The comparison itself is
`_splitwiseinfo(s1) + _splitwiseinfo(s2) - 2 * sum(splitinfo, intersect(s1, s2))`, using
`Splits`'s existing `intersect` (CHUNK-004) and `splitinfo` (CHUNK-013).

After the metric landed, the user asked for it to be benchmarked against TreeDist and added
to `benchmark/`, alongside the existing Robinson-Foulds/quartet/Jaccard-Robinson-Foulds
comparisons. Added `benchmark/inforf.R` (copied from `jrf.R`'s structure: median-of-repeats
timing, GC-peak memory, a `value` column for cross-checking), wired a matching
`benchmarkinforf()` into `run.jl`, and extended `report()` with an Info-Robinson-Foulds
section in `results.md`. A full `julia --project=benchmark benchmark/run.jl` run (with the
Nix R on `PATH`) confirms agreement at every size and surfaces a real finding: this metric
is **10.0×/2.0× faster than TreeDist at n=10/50 but 1.7×/5.6× slower at n=200/1000** —
a steeper falloff than `JaccardRobinsonFoulds`'s (2.6×/2.8× slower) despite doing strictly
less work per pair (no assignment solve, just split intersection and a sum). That is the
signature of the `O(n²)` `log2rooted`/`log2unrooted` recomputation the plan's Open
Questions section already named as a deferred concern — now measured rather than merely
anticipated. Recorded in both `benchmark/README.md` and the plan's Open Questions entry.

The user then asked to restructure the plan around that finding rather than leave it as an
open-ended "revisit later" note. Added **CHUNK-034 (split-information-table)**, depending
only on CHUNK-013, and resequenced the plan's ordering comment (the same informal
"Agreed sequence" idiom used for earlier reprioritizations, not a `Depends on` edit on the
downstream chunks — CHUNK-014/015/020/030/031 aren't functionally blocked by CHUNK-034, they
would just each rediscover the same `O(n²)` problem independently if implemented first) so
CHUNK-034 runs before CHUNK-014/015/020/030/031. The stale Open Questions entry about the
`O(n²)` recomputation is now struck through and points at CHUNK-034 instead of sitting open.

## Key decisions made

- **Citation is Smith (2020), Bioinformatics 36(20):5007–5013, §2.1** — found in TreeDist's
  `inst/REFERENCES.bib` as `SmithDist`, cited from `R/tree_distance_rf.R`'s roxygen block
  for `InfoRobinsonFoulds` specifically (`@insertCite{@§2.1 in @SmithDist}`). Not derived
  from memory alone; read directly from the source this session.
- **A rooted tree matches its unrooted twin exactly (bitwise `0.0`)**, unlike
  `WeightedRobinsonFoulds`. `InfoRobinsonFoulds` only ever touches split *identity*
  (`splitinfo` depends on split size, not branch length), so there is no floating-point
  branch-length summation to introduce a last-ulp disagreement between rootings — this is
  a genuine difference from the branch-length metrics, not an oversight.
- **Cross-check uses a tolerance (`atol=1e-9, rtol=1e-6`), not bitwise equality**, the same
  choice CHUNK-011 made for Nye/JRF and for the same underlying reason stated in the plan's
  Reference implementation section: TreeDist applies `.FloorNumericalNoise` to
  `InfoRobinsonFoulds`'s result, which this package's formula does not attempt to
  reproduce bit-for-bit.

## State of the codebase

- Files modified: `src/robinsonfoulds.jl` (new `InfoRobinsonFoulds` struct, `_compare`,
  `normalizerinfo`, `_splitwiseinfo` helper), `src/PhyloDistances.jl` (new export),
  `test/test_robinsonfoulds.jl` (new `InfoRobinsonFoulds` testset, 56 tests: hand-computed
  value, identical trees, star-vs-resolved, disjoint-split additivity, exact rooted/unrooted
  agreement, branch-length independence, zero-iff-RF-zero and non-negativity over random
  pairs, symmetry, normalization, result type, convention agreement, mismatched-taxa
  rejection), `validation/crosscheck.jl` (added `InfoRobinsonFoulds()` and its normalized
  form to both the R `compare()` function and the Julia comparison, with a tolerance
  comparison and a new `:irf`/`:irfnorm` row in `QUANTITIES`), `ANALYSIS_PLAN.md`
  (CHUNK-033 marked complete with implementation notes, session ledger line, one new Open
  Questions entry — see Watch out for).
- Package loads cleanly: yes, Julia 1.12.6 via the MCP session.
- Test suite passes: yes — `Pkg.test()` reports **7561/7561** (up from 7505 at the start of
  this session, +56 for `InfoRobinsonFoulds`).
- Validated against TreeDist 2.14.1 directly: `validation/crosscheck.jl` (full default
  sweep, 1140 generated cases) reports **0 mismatches** for `InfoRobinsonFoulds()` and
  `InfoRobinsonFoulds(normalize = true)`, alongside the pre-existing zero mismatches for
  every other quantity it checks. `validation/report.md` is updated and staged.
- Additional files this session: `benchmark/inforf.R` (new), `benchmark/run.jl` (new
  `benchmarkinforf()`, `report()` gained a parameter and a section, `main()` wires up the
  new R script), `benchmark/README.md` (new "Info-Robinson-Foulds" section),
  `benchmark/results.md` (regenerated by a full benchmark run), `.gitignore`
  (`benchmark/results_inforf_r.tsv` added to the existing intermediate-output pattern).
- Not yet committed: everything above is staged but not committed, per the project's
  standing preference to show a draft commit message for review first. This is naturally
  two logical commits (the metric itself, then the benchmark addition) if the user wants
  them split.
- Entry point: `InfoRobinsonFoulds()(t1, t2)`, exported. Benchmark entry point:
  `PATH=/nix/store/jaqvbj23b52yl0qgcrrb4ysbxdlqlbv5-R-4.4.2-wrapper/bin:$PATH julia
  --project=benchmark benchmark/run.jl`.
- Known issues: none — the O(n²) slowdown at scale is a known, already-flagged limitation,
  not a bug; see Key decisions made above.

## Next chunk

**CHUNK-034 (split-information-table)**, added and prioritized this session specifically
because of the benchmark finding above. It depends only on CHUNK-013 (complete). Give
`log2rooted`/`log2unrooted` an `O(n)`-total, `O(1)`-per-split table-backed path — computing
`log2rooted(0), ..., log2rooted(n)` once instead of recomputing the running sum per split —
used wherever a metric sums information content over a whole split set (`SplitwiseInfo`).
Keep the existing single-value API unchanged; this adds a second path, not a replacement.
Update `InfoRobinsonFoulds` (CHUNK-033) to use it once it exists, then re-run
`benchmark/inforf.R` to confirm the n=200/1000 regression closes.

After CHUNK-034, the plan's updated sequencing comment (search `ANALYSIS_PLAN.md` for
"Re-prioritized again 2026-08-20") calls for CHUNK-014/015/020/030/031 next, in any order,
since they all build on the same information primitives and should get the fast path from
the start. CHUNK-012 (matching-split-distance) remains independent of all of this — it
depends only on CHUNK-010, not CHUNK-013 — and can be done whenever, unaffected.

## Watch out for

- Everything already flagged in prior handoffs still applies (Nix R path — put
  `/nix/store/jaqvbj23b52yl0qgcrrb4ysbxdlqlbv5-R-4.4.2-wrapper/bin` on `PATH` ahead of
  `/usr/local/bin/R` before running anything in `validation/`; the ~400 uncaptured rooting
  warnings in `test_robinsonfoulds.jl`; `Bool <: Real`; `NewickTree.Node` construction;
  no bare `using AbstractTrees`/`using NewickTree`; give test helpers file-specific names;
  the two hand-verification-repeats-the-bug lessons from CHUNK-011).
- **The scratch TreeDist clone is session-scoped and gone again.** Re-clone with
  `git clone --depth 1 https://github.com/ms609/TreeDist.git` into a fresh scratchpad
  directory next session if source reading is needed — this is now the second session in a
  row this has come up; consider whether a longer-lived clone location is worth it if a
  third session needs the same thing.
- **`validation/README.md` now describes the crosscheck as more exact than it is** — see
  the new Open Questions entry. Its framing prose still says "compared... exactly" /
  "bitwise, no tolerance" while three metric families (Nye/JRF since CHUNK-011,
  InfoRobinsonFoulds since this session) are actually compared to a tolerance. Cheap fix,
  not done this session per scope discipline.
- `InfoRobinsonFoulds`, like `RobinsonFoulds`, excludes trivial splits by default
  (`splits(tree, index)` with no `trivial = true`), which is harmless for the information
  sum specifically — a trivial split's `splitinfo` is exactly zero regardless — but keeps
  the split *count* semantics consistent with classic RF, which is what makes the
  star-vs-resolved-tree test ("costs the resolved tree's full `SplitwiseInfo`") come out
  exact rather than needing to account for pendant-branch information that happens to be
  zero anyway.
