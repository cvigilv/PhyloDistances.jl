# Session Handoff — 2026-08-17

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-002: design-distance-api

`src/interface.jl` defines the surface every metric must fit, built on the **Distances.jl**
interface rather than a private one: `TreeMetric <: Distances.SemiMetric`, applied by
calling it, `m(t1, t2)`. `pairwise`, `pairwise!`, `colwise`, `colwise!`, `evaluate` and
`result_type` all work unchanged and are re-exported. Alongside it are a `Convention`
hierarchy selecting between published formulations, the `requiresrooted` trait, and
`TreeSimilarity` for quantities that are largest on identical trees. No concrete metric
exists yet — this chunk is the contract, verified with dummy metrics defined in the tests.

**Adding a metric requires:** `PhyloDistances._compare(m, ::TreeDistConvention, t1, t2)`,
`convention` and `normalize` fields (or overrides of `convention`/`isnormalized`), and
`Distances.result_type` if the result is not `Float64`. Optional: a `PrimaryConvention`
method, `normalizer`, `requiresrooted`.

## Key decisions made

- **Distances.jl is the interface, not a model to copy.** An earlier iteration reimplemented
  `PreMetric`/`SemiMetric`/`pairwise` badly; that was replaced by the real dependency.
  Metrics are applied by calling them and there are **no call keywords** — this is forced,
  not stylistic, because `pairwise` applies a metric without keywords, so anything not on
  the type is unreachable through it.
- **Every variant selector is a field**: `convention`, `normalize`, and metric-specific
  parameters such as JRF's `k` or Kendall-Colijn's λ. A constructed metric specifies the
  computation completely.
- **Similarities sit outside the Distances hierarchy**, subtyping `TreeSimilarity` with
  their own `Distances.pairwise` methods. `PreMetric` requires `d(x,x) == 0`, and
  `SemiMetric` is actively dangerous: `pairwise` *assumes* the zero diagonal rather than
  computing it. Verified — a similarity declared `SemiMetric` yields
  `[0 98 95; 98 0 97; 95 97 0]` where the truth is `[100 98 95; 98 100 97; 95 97 100]`.
- **`Distances.result_type` must be defined by every concrete metric** whose result is not
  `Float64`. Its default infers from `one(eltype(a))`, which assumes numeric observations
  and throws outright for trees.
- **Distances' functions are re-exported**, so `using PhyloDistances, Distances` is clean:
  the names refer to identical bindings, which cannot clash. Verified.
- **`Distances.pairwise` discards offset axes** — results are 1-based and index by position.
  Accepted rather than diverging from upstream; a test pins the behavior so a change
  upstream is noticed.
- **`TreeMetric <: Distances.SemiMetric`, not `Metric`** — the safe common denominator,
  since not every tree distance satisfies the triangle inequality. Promoting individual
  metrics is an open question.
- **No rooting stub.** Only the `requiresrooted` trait is defined here. CHUNK-003 adds the
  actual precondition check; there is deliberately no placeholder to mistake for working
  code.

## State of the codebase

- Files created: `src/interface.jl`, `test/test_interface.jl`
- Files modified: `src/PhyloDistances.jl`, `test/runtests.jl`, `test/Project.toml`,
  `Project.toml`, `README.md`
- Dependencies added: `Distances` (main), `Distances` + `OffsetArrays` (test-only)
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — 59 tests
- Entry points: no analysis entry point yet.
  `julia --project=. -e 'using Pkg; Pkg.test()'` runs the suite.
- Known issues: none

## Next chunk

CHUNK-003: tree-ingest-and-taxa

Newick ingest via NewickTree.jl plus canonical taxon indexing: the sorted leaf-label vector
and a label→index map for one tree, and a check that two trees span the same taxon set that
throws naming the differing labels. This chunk also implements the rooting policy: a
predicate detecting whether a parsed tree is rooted, and a precondition helper that compares
it against `requiresrooted(metric)`, **warns** naming tree, metric and conversion, and
applies that conversion — or throws where no defensible conversion exists.

Wire that helper into `_apply` in `src/interface.jl`, which is the single funnel every
metric passes through.

All traversal must go through AbstractTrees.jl, not `NewickTree.Node` internals, so the
package accepts any conforming node type.

## Watch out for

- **The specific rooting conversion rules are still an open question**, and CHUNK-003 is
  where they get settled. Unrooted → rooted needs a documented rooting rule (midpoint?
  outgroup-required?); rooted → unrooted needs a root-split policy. Whatever is chosen is
  visible in every warning the package emits, so record it in Working knowledge.
- **The warn-and-convert policy is a deliberate exception to the project's fail-fast
  stance.** It is only defensible because the warning states exactly what conversion was
  applied. Where no such statement can be made, throw.
- **The critical path is CHUNK-001 → CHUNK-008**, ending in validated Robinson-Foulds and
  quartet distance, which are needed for immediate use. Nothing past CHUNK-008 is urgent.
- **Dummy metrics in `test/test_interface.jl` use integers as stand-in trees**, which is
  sound only because the interface performs no tree operations. Once CHUNK-003 wires the
  rooting check into `_apply`, those dummies will need real trees or a rooting predicate
  that tolerates them — expect to revisit that file.
- **`test/fixtures/` and `scripts/` hold `.gitkeep` placeholders.** Delete each when the
  directory gains real content.
- **Julia 1.12 is the declared floor and there is no LTS back-support.** Use 1.12 features
  freely; do not add 1.10 compatibility shims.
