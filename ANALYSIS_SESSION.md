# Session Handoff — 2026-08-17

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-001: scaffold-package

The repository was turned into a working Julia package. PkgTemplates.jl generated the
skeleton in a scratch directory (it refuses to write into an existing directory), and the
wanted files were copied in so the pre-existing `LICENSE`, `README.md`, and `.gitignore`
survived untouched. NewickTree and AbstractTrees were added as dependencies with compat
bounds, and the package now loads and tests cleanly on Julia 1.12.6.

## Key decisions made

- **Julia 1.12 is the minimum supported version; the 1.10 LTS is not supported.** Language
  features and stdlib methods up to 1.12 may be used freely, and `julia +lts` is not part
  of the verification loop.
- **Separate test environment.** Tests resolve against `test/Project.toml` rather than the
  `[extras]`/`[targets]` mechanism. Test-only dependencies (StableRNGs, CairoMakie, …) go
  there, not in the top-level `Project.toml`.
- **Disabled generator plugins.** License, Readme, CompatHelper, TagBot, and Dependabot were
  turned off — the first two would have overwritten existing files, the last three are
  registration concerns that do not apply to an unregistered package.
- **CI matrix is `1.12`, `1`, `pre`** — the declared floor, the current stable release, and
  prereleases.
- **Compat bounds are the newest registered versions** (`AbstractTrees = "0.4.5"`,
  `NewickTree = "0.4.0"`) because both packages are at their latest release. Revisit if a
  lower bound ever needs to be relaxed for a downstream consumer.

## State of the codebase

- Files created: `Project.toml`, `src/PhyloDistances.jl`, `test/runtests.jl`,
  `test/Project.toml`, `test/fixtures/.gitkeep`, `scripts/.gitkeep`,
  `.github/workflows/CI.yml`
- Files preserved unchanged: `LICENSE`, `README.md`, `.gitignore`
- Package loads cleanly: yes, on 1.12.6
- Test suite passes: yes (the suite is empty — this is a scaffolding chunk)
- Entry points: none yet. `julia --project=. -e 'using Pkg; Pkg.test()'` runs the suite.
- Known issues: none

## Next chunk

CHUNK-002: design-distance-api

Establish the public API before any metric exists. It must build in two things that the
settled design decisions require and that would be painful to retrofit across fifteen
metrics: a `convention` keyword (`:treedist` default, `:primary` opt-in) threaded through
dispatch, and an `isrooted(::Metric)` trait that CHUNK-003's precondition check consults to
decide whether to warn and convert. Also ship a naive all-pairs fallback (loop over pairs)
so a distance matrix works from day one; CHUNK-024 optimizes it without changing the API.

Verify with a dummy metric defined inside the test suite, exercised through every documented
entry point and both `convention` values.

## Watch out for

- **The critical path is CHUNK-001 → CHUNK-008.** RF and quartet distance are needed for
  immediate use. CHUNK-002 sets the API surface all fifteen metrics must fit, so it is worth
  getting right, but do not gold-plate it — nothing past CHUNK-008 is urgent.
- **`test/fixtures/` and `scripts/` hold `.gitkeep` placeholders.** Delete each when the
  directory gains real content.
- **`src/PhyloDistances.jl` has `using AbstractTrees` and `using NewickTree` but exports
  nothing.** Decide the export policy in CHUNK-002 rather than accumulating exports
  incidentally.
- **Julia 1.12 is the declared floor and there is no LTS back-support.** Use 1.12 features
  freely; do not add 1.10 compatibility shims.
