# Session Handoff — 2026-08-19

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-018: quartet-distance-fast. `QuartetDistance` now defaults to an `O(n³)` algorithm
that counts concordant quartets without enumerating them, replacing the `O(n⁴)` direct
enumeration as the default (that enumeration is still available via `algorithm = :naive`,
and remains the correctness oracle the fast path is tested against).

At n=1500 (a size the user specifically wanted usable), a single pairwise call dropped from
an extrapolated ~15 minutes to ~7–8 seconds — checked on both a random binary tree and a
100-taxon caterpillar (the shape most likely to expose an unbalanced-tree bug).

## Key decisions made

- **O(n³), not tqDist's O(n log n).** Asked the user directly: tqDist's algorithm needs a
  dynamic Hierarchical Decomposition Tree ported from a research paper, high implementation
  risk with a failure mode (silent wrong answers on large unbalanced trees) that an n≤12
  brute-force cross-check would not reliably catch. The user chose the lower-risk O(n³)
  trade explicitly. If ~7 s per pair is later found inadequate, tqDist is the documented
  next step (Open Questions).
- **The fast path requires at least one input tree to be binary.** For two polytomous
  trees, exactness would need a materially harder computation (see Working knowledge in the
  plan). Rather than implement that or silently give a wrong answer, `:fast` warns and falls
  back to `:naive` — a real performance cliff for that case, but a visible one, consistent
  with the project's fail-fast stance. Confirmed both directions work: the fast path handles
  one polytomous tree exactly (no fallback needed, since an unresolved quartet there can
  never be concordant against a binary tree), and both-polytomous correctly triggers the
  warning and matches `:naive`.
- **`algorithm` is a field, not a call keyword**, per the project's settled design decision;
  `QuartetDistance(; algorithm = :fast)` (default) or `:naive`, validated at construction.

## State of the codebase

- Files modified: `src/quartet.jl` (the fast algorithm, `_naivequartetdistance` split out of
  the old `_compare` body, the `algorithm` field), `src/clustertable.jl` (`_rootedorder`
  gained an optional `target` taxon-position parameter, default `1`, so it can root a
  `FlatTree` at any leaf — Robinson-Foulds' call site is unchanged), `test/test_quartet.jl`
  (algorithm-selection and fast-vs-naive tests).
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — **1319 tests**, up from 1274, no R and no network
- Entry points: `QuartetDistance()(t1, t2)` (fast by default),
  `QuartetDistance(; algorithm = :naive)(t1, t2)` (the O(n⁴) oracle)
- Known issues: none in what's implemented; the two-polytomy case is a documented
  performance gap, not a bug (see Open Questions in the plan)

## Algorithm, briefly

An unrooted quartet's topology is a rooted triple under any one of its four members as
outgroup, and every quartet is discovered exactly four times this way (once per member).
Summing rooted-triple agreement between the two trees over all `n` choices of root and
dividing by four gives the quartet agreement count — turning
"enumerate `binomial(n,4)` quartets" into "for each of `n` roots, count agreement over
`binomial(n-1,3)` triples without enumerating those either." For a fixed root, leaf pairs
are grouped by their most-recent-common-ancestor in tree 1 (discovered for free while
walking tree 1's structure), and for each group a T2-position-indexed membership array
(built once per branch point, `O(n)`) turns "how many of this T1 clade's leaves lie in a
given T2 clade" into an `O(1)` prefix-sum lookup — giving `O(n²)` per root, `O(n³)` overall.
Full derivation and the two-position-numbering trap that cost an hour of debugging are in
the plan's Working knowledge, dated 2026-08-19 under CHUNK-018.

## Verification beyond the plan's bar

- Naive-vs-fast cross-check: 270 random binary-tree pairs (n=4–12), 30 more at n=20/50/100,
  240 one-polytomy pairs and 240 both-polytomy pairs (n=4–15), a 100-taxon caterpillar pair.
  All exact matches. The n=1500 and caterpillar timing runs were exploratory (not committed);
  the committed suite (`test/test_quartet.jl`) keeps a representative slice.
- An internal `total % 4 == 0` assertion inside `_fastconcordantcount` (every concordant
  triple must be discovered at exactly 4 of the `n` roots) caught the position-numbering bug
  described above before it reached the naive-comparison tests — worth keeping as a
  standing correctness check, not just a debugging aid.

## Next chunk

Per the agreed sequence (Working stance note, set 2026-08-19): **CHUNK-010
(split-matching-framework)** next — the vendored optimal-assignment solver and the
generalized-RF matching engine that CHUNK-011 through CHUNK-015 reduce to. CHUNK-013
(split-information-primitives) is independent and can be built before or alongside it.

## Watch out for

- Everything already flagged in prior handoffs still applies (Nix R path, `readnw` and
  `SubString`, the ~400 uncaptured rooting warnings in `test_robinsonfoulds.jl`, `Bool <:
  Real`, `NewickTree.Node` construction, no bare `using AbstractTrees`/`using NewickTree`).
- **This worktree's `ANALYSIS_PLAN.md` was rebased from the shared checkout's *last commit*
  (085c99d), not its working tree.** The shared checkout had uncommitted edits to the plan
  (the "Agreed sequence" note, richer CHUNK-018/CHUNK-010 notes) that were never committed;
  this session read them from the shared checkout directly and folded them into this
  worktree's copy alongside the CHUNK-018 completion notes, so nothing should be lost — but
  when merging, diff against the shared checkout's current `ANALYSIS_PLAN.md` rather than
  assuming a plain fast-forward captures everything.
- If a future session extends the fast algorithm to handle two polytomous trees, the
  building block it is missing is a per-node-pair contingency-table count (leaves split
  across ≥3 children of a branch point in *both* trees simultaneously) — sketched but not
  derived in detail in the plan's Open Questions.
