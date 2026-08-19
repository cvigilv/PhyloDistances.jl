# Session Handoff — 2026-08-19

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-010 (split-matching-framework) and CHUNK-011 (nye-similarity-and-jrf), done together
in one session at the user's request ("implement JRF and validate against TreeDist").

CHUNK-010 built the generalized Robinson-Foulds reduction: a vendored linear-assignment
solver (`src/assignment.jl`, `_hungarian` — the classical shortest-augmenting-path
algorithm with vertex potentials, `O(min(n,m)² max(n,m))`, working on exact `Float64`
costs) and a public matching framework (`src/splitmatching.jl`, `splitmatching(scorer,
splits1, splits2; maximize)`) that turns any pairwise split-scoring function into an
optimal-matching total.

CHUNK-011 built `NyeSimilarity` (a `TreeSimilarity`) and `JaccardRobinsonFoulds` (a
`TreeMetric`, fields `k` and `allowconflict`) on top of that framework, both scored by
`_jaccardscore` — a port of TreeDist's `jaccard_similarity` C++ scorer
(`src/generalizedrf.jl`).

## Key decisions made

- **Both Nye and JRF, not JRF alone.** Asked the user; they chose "both" since Nye is JRF's
  `k=1, allowconflict=true` case computed as a similarity, so building JRF builds Nye's
  scoring for free.
- **The vendored solver is Kuhn-Munkres/Munkres over exact `Float64`, not a port of
  TreeDist's integer-quantized Jonker-Volgenant.** TreeDist scales costs to `Int64` (`BIG =
  typemax(Int64) ÷ SL_MAX_SPLITS`) for its solver; this package doesn't need that trick in
  Julia and uses exact floating-point costs throughout. Consequence: values agree with
  TreeDist to a tolerance (`atol=1e-9, rtol=1e-6`), not bitwise — see Working knowledge.
- **`splitmatching` is `public` but not exported**, matching the treatment of
  `normalizerinfo`, `convention`, etc. — an extension point, not an everyday call.

## State of the codebase

- Files created: `src/assignment.jl`, `src/splitmatching.jl`, `src/generalizedrf.jl`,
  `test/test_assignment.jl`, `test/test_splitmatching.jl`, `test/test_generalizedrf.jl`,
  `benchmark/jrf.R`.
- Files modified: `src/PhyloDistances.jl` (new includes/exports), `test/runtests.jl` (new
  includes), `validation/crosscheck.jl` (Nye/JRF cross-check against TreeDist, extended
  `checkpair`/`QUANTITIES`/report text), `benchmark/run.jl` (`benchmarkjrf`, `jrf.R`
  invocation, new report section, banner deduplication), `.gitignore`
  (`benchmark/results_jrf_r.tsv`).
- Package loads cleanly: yes, Julia 1.12.6.
- Test suite passes: yes — **4448 tests**, up from 4444 (net +4 after adding several
  hundred new tests and removing nothing), no R and no network for the committed suite.
- Entry points: `NyeSimilarity()(t1, t2)`, `JaccardRobinsonFoulds()(t1, t2)`,
  `JaccardRobinsonFoulds(; k = 2, allowconflict = false)(t1, t2)`,
  `PhyloDistances.splitmatching(scorer, splits1, splits2)` for a custom metric.
- Known issues: none functionally; a real performance gap is documented (below).

## Validation against TreeDist

`validation/crosscheck.jl` (`julia --project=validation validation/crosscheck.jl 3000`,
run this session with Nix R 4.4.2 on `PATH`): **3,140 generated cases, 0 mismatches** for
`NyeSimilarity()`, `NyeSimilarity(normalize=true)`, `JaccardRobinsonFoulds()`,
`JaccardRobinsonFoulds(k=2, allowconflict=false)`, and
`JaccardRobinsonFoulds(normalize=true)`, alongside the pre-existing RF/quartet checks
(also 0 mismatches). Committed report: `validation/report.md`.

Unlike RF and the quartet distance, this comparison is **tolerance-based, not bitwise** —
see Working knowledge for why that is the correct bar here, not a weaker one.

## Benchmark against TreeDist

`benchmark/run.jl` (run this session): **PhyloDistances is currently slower than TreeDist
for this metric, and the gap widens with tree size** — 4.5× faster at n=10, but 2.0×
slower at n=50, 5.9× slower at n=200, 7.0× slower at n=1000. This is the opposite of the
RF and quartet-distance results, both of which are faster here at every measured size.
Diagnosed (not yet fixed — out of this chunk's scope) in `benchmark/README.md`: the
`Splits` type stores each split as a `BitVector`, so `splitmatching`'s O(n²) score-matrix
build does one allocating `count(a .& b)` per cell, where TreeDist packs splits into
machine words and computes the same overlap with masked integer ops. Full numbers in
`benchmark/results.md`.

## Next chunk

Per the plan's session ledger: **CHUNK-012 (matching-split-distance)**, which depends only
on CHUNK-010 (now complete) — a straightforward second scoring function over
`splitmatching`. CHUNK-013 (split-information-primitives) remains independent and can be
picked up instead or alongside.

**Before either**, consider whether to spend a chunk on the `Splits` `BitVector` → packed
representation switch flagged above — CHUNK-012, CHUNK-014, and CHUNK-015 all build more
split-matching metrics on the same representation, so the performance gap will replicate
across each of them if left unaddressed. Not urgent (correctness is unaffected either
way), but cheaper to fix once than to duplicate.

## Watch out for

- Everything already flagged in prior handoffs still applies (Nix R path, `readnw` and
  `SubString`, the ~400 uncaptured rooting warnings in `test_robinsonfoulds.jl`, `Bool <:
  Real`, `NewickTree.Node` construction, no bare `using AbstractTrees`/`using NewickTree`).
- **Test files share one namespace.** `test/runtests.jl` `include`s every file into one
  `@testset` in `Main`, not separate modules. A same-name helper in two files silently
  shadows rather than errors, and if the arities coincide too (`caterpillar(x)` in two
  files, taking different kinds of `x`), the second file's definition permanently replaces
  the first's for anything included afterward. This session renamed collisions found this
  way (`mask` → `splitmask`, `caterpillar` → `orderedcaterpillar` in the new files) but did
  not audit the rest of the suite for the same risk — worth a pass if it ever causes a
  confusing failure.
- **A hand-verification that repeats the same arithmetic the code performs will not catch
  a bug in that arithmetic.** This session's `_jaccardscore` bug (`nA` instead of `nB` in
  one term) passed hand-checked cases for exactly this reason — the check used the same
  wrong formula. It was caught by matching-order symmetry, a property independent of how
  the scorer computes its answer. Prefer structural checks (symmetry, a known limit,
  brute-force) over redoing an implementation's derivation by hand when verifying a port.
