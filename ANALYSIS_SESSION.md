# Session Handoff — 2026-08-17

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-002: design-distance-api

`src/interface.jl` defines the surface every metric must fit: a `TreeMetric` supertype, a
`Convention` hierarchy selecting between published formulations, the `requiresrooted` and
`issimilarity` traits, `compare`, and a naive `pairwise`. No concrete metric exists yet —
this chunk is the contract, verified with dummy metrics defined inside the test suite.

**Adding a metric requires exactly one method:**
`PhyloDistances._compare(m, ::TreeDistConvention, t1, t2)`. Every other part of the
interface has a fallback, and the two that can legitimately be missing (`:primary`
convention, `normalizer`) fail with a message naming the metric rather than a `MethodError`.

## Key decisions made

- **`compare(metric, t1, t2)` is the single entry point.** Metrics are deliberately *not*
  callable — `m(t1, t2)` does not work, and a test asserts no call method exists. One way to
  apply a metric, not two.
- **Exported: `TreeMetric`, `Convention`, `TreeDistConvention`, `PrimaryConvention`,
  `compare`, `pairwise`.** The traits (`normalizer`, `requiresrooted`, `issimilarity`) are
  `public` but unexported (Julia 1.11+ `public` keyword) since they are extension points
  rather than everyday calls. `pairwise` collides with `Distances.pairwise` and
  `StatsBase.pairwise`; it is exported anyway and the collision is deferred — see the plan's
  Open Questions. Changing this after release would be breaking.
- **Metric parameters are fields of the metric type**, not call keywords — so JRF's `k` and
  Kendall-Colijn's λ live on the type, and a constructed metric specifies the computation
  completely. `convention` and `normalize` are the only keywords, which is what lets
  `pairwise` forward everything without knowing about any specific metric.
- **`convention` accepts a symbol or a `Convention` instance.** `Convention(::Symbol)`
  resolves `:treedist`/`:primary` at the public boundary so internals dispatch on singleton
  types; `Base.Symbol(::Convention)` gives the reverse for error messages.
- **Normalization is a division.** A metric defines `normalizer(m, conv, t1, t2)` and
  `compare` divides by it. A metric with no normalization throws rather than silently
  returning the raw value.
- **`pairwise` honors generic axes** — the result's axes track the input's, tested against
  `OffsetVector`. It fills one triangle and mirrors, relying on the documented symmetry of
  every metric. Element type is taken from the diagonal, which is part of the result, so no
  comparison is computed and discarded.
- **`OffsetArrays` added as a test-only dependency** (`test/Project.toml`) to enforce the
  generic-indexing contract.
- **No rooting stub.** Only the `requiresrooted` trait is defined here. CHUNK-003 adds the
  actual precondition check into `compare`; there is deliberately no placeholder to
  mistake for working code.

## State of the codebase

- Files created: `src/interface.jl`, `test/test_interface.jl`
- Files modified: `src/PhyloDistances.jl`, `test/runtests.jl`, `test/Project.toml`,
  `README.md`
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — 54 tests
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
applies that conversion — or throws where no defensible conversion exists. Wire that helper
into `compare` in `src/interface.jl`.

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
- **`test/fixtures/` and `scripts/` hold `.gitkeep` placeholders.** Delete each when the
  directory gains real content.
- **Julia 1.12 is the declared floor and there is no LTS back-support.** Use 1.12 features
  freely; do not add 1.10 compatibility shims.
- **Dummy metrics in `test/test_interface.jl` use integers as stand-in trees**, which is
  sound only because the interface performs no tree operations. Once CHUNK-003 wires the
  rooting check into `compare`, those dummies will need real trees or a metric declaring
  `requiresrooted` in a way the check tolerates — expect to revisit that file.
