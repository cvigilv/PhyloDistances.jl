# Session Handoff — 2026-08-17

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-007: quartet-distance-exact, plus the reference validation and benchmarks for it.

`src/quartet.jl` adds `QuartetDistance`: the number of four-taxon subsets the two trees
resolve differently, by direct enumeration of all `binomial(n, 4)` quartets. A quartet's
topology is read off **leaf-to-leaf path lengths** rather than off the split set.
`_topologicalpaths` tabulates the edge count between every pair of leaves in O(n²) — each
pair filled once, at the node where its two subtrees join — and each quartet is then resolved
in constant time by the four-point condition: of the three sums of opposite-pair distances,
the smallest names the topology, and a three-way tie means the tree resolves it not at all.

The R package **Quartet 1.3.0 was installed** to serve as its reference, and the validation
and benchmark harnesses were extended to cover it.

## Key decisions made

- **A quartet resolved by only one tree counts as a difference** — confirmed to match the
  reference exactly, not merely chosen. Quartet's `QuartetStatus` reports counts rather than
  a distance, and ours is its `d + r1 + r2`.
- **`normalize = true` divides by `binomial(n, 4)`**, which the crosscheck holds equal to
  Quartet's own `Q`, so both sides agree on the divisor as well as the distance.
- **Path lengths, not splits.** Resolving a quartet by scanning the split set would cost a
  factor of n more and is no simpler.
- **Deliberately the naive O(n⁴) algorithm**, as the chunk specifies. It is slower than the
  reference and that is expected: tqDist counts the same quantity without enumerating it.
- **R is reached two different ways, on measured grounds** — see below.

## The harness split: RCall for validation, subprocess for benchmarks

This was evaluated with measurements rather than taste, and the reasoning should survive.

**Validation uses RCall** (`validation/crosscheck.jl`). Values cross as machine doubles and
integers, so nothing is formatted. That *retires* the old "compare in the direction R →
Julia" rule rather than following it: that rule existed only because R's `as.numeric` cannot
reliably round-trip its own `%.17g` output. With no text there is no round trip to get wrong.
R's `stop()` also becomes a catchable Julia exception and `NA_integer_` becomes `missing`.

**Benchmarks keep R in a subprocess** (`benchmark/*.R`). An RCall round trip costs **12 µs**
bare and **30 µs** with an argument, against about **1 µs** for the same call timed inside R.
TreeDist's Robinson-Foulds on ten taxa is 37 µs in total, so the overhead would be a third to
four fifths of the measurement. The timing loop must live inside R either way, which leaves
RCall adding a dependency and a failure mode for nothing. A subprocess also gives each R
benchmark a cold interpreter, keeping the `gc()` figures honest.

`ENV["R_HOME"] = strip(read(\`R RHOME\`, String))` must be set **before** `using RCall`: RCall
records R's location when built, and a Nix rebuild moves that store path.

## State of the codebase

- Files created: `src/quartet.jl`, `test/test_quartet.jl`, `benchmark/quartet.R`
- Files modified: `src/PhyloDistances.jl`, `test/runtests.jl`, `.gitignore`,
  `benchmark/run.jl`, `benchmark/README.md`, `benchmark/results.md`,
  `validation/crosscheck.jl`, `validation/README.md`, `validation/report.md`,
  `validation/Project.toml`
- Files deleted: `validation/treedist_values.R` (the crosscheck no longer shells out)
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — 963 tests, no R required
- Entry points: `QuartetDistance()(t1, t2)`, `pairwise(QuartetDistance(), trees)`
- Known issues: none

## Usable right now

```julia
using PhyloDistances
t1 = readnw("(((Human,Chimp),Gorilla),Orang,Gibbon);")
t2 = readnw("(((Human,Gorilla),Chimp),Gibbon,Orang);")

RobinsonFoulds()(t1, t2)                      # split-symmetric-difference distance
QuartetDistance()(t1, t2)                     # four-taxon subsets resolved differently
QuartetDistance(; normalize = true)(t1, t2)   # as a fraction of all C(n,4) quartets
pairwise(QuartetDistance(), [t1, t2, ...])    # all-pairs matrix
```

## Validation status

`julia --project=validation validation/crosscheck.jl` — **740 cases, zero mismatches** across
five quantities: `RobinsonFoulds` raw and normalized against TreeDist 2.14.1, and
`QuartetDistance` raw and normalized plus the quartet count `Q` against Quartet 1.3.0.
Integers exactly, floats bitwise, `NaN` meeting `NaN`. It exits non-zero on any difference,
so it can gate a release. `validation/report.md` is the committed record.

The test suite stays portable and R-free: quartet distance is covered there by hand-computed
cases and an **independent split-based oracle** (a tree resolves `ab|cd` iff some split
separates that pair).

## Benchmark status

`julia --project=benchmark benchmark/run.jl` — R stays in a subprocess, both sides read the
same Newick files, neither timing includes parsing. `benchmark/results.md` is the committed
record and `run.jl` asserts both implementations returned the same distance for every pair it
timed.

| taxa | quartets | PhyloDistances | Quartet | ratio |
|-----:|---------:|---------------:|--------:|------:|
| 10 | 210 | 9.7 µs | 1.71 ms | 176.4× faster |
| 50 | 230,300 | 821.2 µs | 12.21 ms | 14.9× faster |
| 200 | 64,684,950 | 222.34 ms | 183.46 ms | 1.2× slower |
| 477 | 2,130,031,575 | 7.68 s | 935.14 ms | 8.2× slower |
| 700 | 9,918,641,075 | 37.40 s | — | — |
| 1000 | 41,417,124,750 | 183.30 s | — | — |

**The crossover is near 200 taxa.** Below it R's per-call overhead dominates and this package
wins comfortably; above it the algorithms decide — exact `O(n⁴)` enumeration here against
tqDist, which counts without enumerating — and the gap widens without bound. That is the
measured case for CHUNK-018, and the reason the docstring states the complexity plainly.

Robinson-Foulds is unchanged and still 5–29× faster than TreeDist on a single pair, 3.1×
slower on all-pairs (the outstanding CHUNK-024 work).

## Traps found this session, all still live

- **TreeDist and Quartet both export `RobinsonFoulds`**, meaning different things. Whichever
  is attached second wins, and Quartet's is a deprecated function expecting a status vector,
  so it fails deep inside with `subscript out of bounds`. Namespace-qualify every call in R
  code that touches both: `TreeDist::RobinsonFoulds`, `Quartet::QuartetStatus`.
- **Quartet's 477-tip ceiling bites just below the refusal.** Above 477 it stops with "trees
  too large for integer representation"; *at* 477 the `N` column is already `2 * Q` overflowed
  to `NA` while `Q` itself still fits. Read `Q`/`s`/`d`/`r1`/`r2`/`u`, never `N`.
- **`QuartetStatus` returns counts, not a distance.** Deciding which counts to combine *is*
  choosing the convention.

## Next chunk

CHUNK-008: validate-rf-and-quartet — the "usable today" checkpoint that closes the critical
path.

Much of its original purpose is now met by `validation/`, so its remaining job is the
**portable, committed fixture**: seed `test/fixtures/` with tree pairs (identical, one NNI
apart, maximally different, a rooted input exercising the warning path, a polytomy,
zero-length branches), each with both metrics' expected values and a provenance line, and the
tests that check against it. Values can now be sourced from the R references with confidence,
but the fixture itself must be plain text with no absolute paths, and the tests must pass
with no R and no network. Also confirm the naive all-pairs API gives a sensible matrix over a
small collection.

## Watch out for

- **`test/fixtures/` still holds a `.gitkeep` placeholder.** CHUNK-008 gives it real content;
  delete the placeholder then. `scripts/` likewise, later.
- **A full benchmark run takes several minutes**, nearly all of it the quartet distance at
  700 and 1000 taxa. Do not run it alongside anything else, or the timings are meaningless.
- **A node with one branch below it inflates path lengths across it** — a rooted tree's root
  is one. Harmless for quartets (subdivision preserves the tree metric), but *not* harmless
  for anything counting nodes or clusters: it silently inflated every RF distance involving a
  rooted input in CHUNK-006. Treat each new traversal on its own terms.
- **Randomized checks must include rooted and polytomous trees.** `randomtree` produces
  neither; 2,000 random unrooted pairs missed the CHUNK-006 root bug that seven hand-picked
  pairs caught. Polytomies are also the only place the two quartet conventions differ.
- **Check whether a reference implementation exists before assuming one.** TreeDist has no
  weighted RF and no quartet distance; the branch-length metrics still have no reference at
  all until phangorn is installed (CHUNK-032).
- **Derive normalizers from the reference, never from the literature** — TreeDist normalizes
  RF by `n1 + n2`, not `2(n-3)`.
- **`.FloorNumericalNoise` zeroes TreeDist's values** below `sqrt(eps) * max(1, info)`.
  Metrics here implement the plain formula; measure the divergence during validation.
- **`Bool <: Real` in Julia**, so `normalize = true` means the metric's own scheme and
  `normalize = 1` means divide by one.
- **Length arithmetic across rootings is inexact**, so branch-length metrics assert with `≈`.
- **`splits` excludes trivial splits by default**; pass `trivial = true` to sum over all
  branches.
- **`NewickTree.Node(id, data)` leaves `parent` and `children` undefined**; build with
  `push!(parent, child)`.
- **Do not reintroduce a bare `using AbstractTrees` or `using NewickTree`** — both export
  `children`, `getroot`, `isroot` and `print_tree`.
- **Julia 1.12 is the declared floor; there is no LTS back-support.** No Julia MCP server was
  available this session, so everything ran as one-shot scripts through Bash.
