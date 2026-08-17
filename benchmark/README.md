# Benchmarks

Compares this package against the R packages whose values it reproduces —
[TreeDist](https://github.com/ms609/TreeDist) for Robinson-Foulds and
[Quartet](https://github.com/ms609/Quartet) for the quartet distance — so that agreeing on
results does not quietly cost an order of magnitude in speed.

## Running

```console
$ julia --project=benchmark benchmark/run.jl
```

`run.jl` writes Newick files to `benchmark/trees/`, benchmarks this package, invokes
`treedist.R` and `quartet.R` on the same files, and renders `results.md`. R is optional —
without a working `Rscript` the Julia timings are reported alone.

**A full run takes several minutes**, nearly all of it the quartet distance at 700 and 1000
taxa, where a single call is a minute or two.

Neither side's timing includes parsing. `benchmark/trees/` and the `results_*_r.tsv` files
are regenerated on every run and are not tracked; `results.md` is.

Memory is measured differently on each side and the two figures are **not** comparable.
Julia reports bytes allocated by the call; R reports the peak its garbage collector saw,
which includes everything already resident.

## Why R runs as a subprocess and not through RCall

The validation harness calls R through [RCall](https://github.com/JuliaInterop/RCall.jl),
which is the better tool there. It is the wrong tool here.

An RCall round trip costs **12 µs** bare and **30 µs** with an argument, against roughly
**1 µs** for the same call timed inside R. TreeDist's Robinson-Foulds on ten taxa takes 37 µs
in total, so the harness overhead would be a third to four fifths of the measurement. The
timing loop has to run inside R whichever way R is invoked, which leaves RCall adding a
dependency and a failure mode without buying anything.

A subprocess also gives each R benchmark a cold interpreter, which keeps the `gc()` figures
from accumulating whatever a long-lived embedded R happened to be holding.

## Reading the results

The picture at the time of writing, from `results.md`.

### Robinson-Foulds

| taxa | PhyloDistances | TreeDist | ratio |
|-----:|---------------:|---------:|------:|
| 10 | 1.5 µs | 43.9 µs | 28.8× faster |
| 50 | 5.4 µs | 53.2 µs | 9.9× faster |
| 200 | 24.5 µs | 129.0 µs | 5.3× faster |
| 1000 | 111.2 µs | 2.56 ms | 23.0× faster |
| all pairs, 40 trees × 60 taxa | 6.13 ms | 2.00 ms | 3.1× slower |

Robinson-Foulds is computed by Day's (1985) cluster-table algorithm, the same one TreeDist
uses. Rooting a tree at one taxon and numbering its leaves depth-first makes every cluster a
contiguous run of numbers, so a cluster is a pair of endpoints and membership is a
constant-time lookup in a table with two rows per taxon. Comparing two trees then costs time
proportional to their size, where comparing split sets costs the product of their split
counts.

Two further things mattered as much as the algorithm, and neither was visible without
profiling:

**Reading a tree once instead of three times.** Deriving the taxon ordering, validating it
and then reading the structure each walked the tree separately, and iterating
`AbstractTrees.Leaves` allocates a cursor per step — 39,700 of them for a single comparison
at 200 taxa. One typed, iterative walk now collects structure and labels together, taking
`flatten` from 45.7 µs to 3.5 µs and a whole comparison from 5,346 allocations to 172.

**Not sorting labels the answer does not depend on.** Comparing clusters needs only that
both trees agree on which taxon is which, never that the numbering means anything outside
the comparison. Sorting — which a `TaxonIndex` does so that split masks stay reproducible —
was 64% of the remaining time at 1000 taxa, spent hashing and comparing strings. Numbering
the second tree against the first instead halved the total.

**`@inbounds` is deliberately absent.** Profiling put the cost in string hashing and
allocation, never in array indexing: the encoding's own loops are a minority of the runtime
and only one frame of 1,399 showed runtime dispatch. Removing bounds checks would trade
silent undefined behaviour for a few percent of something that is not the bottleneck.

**The all-pairs gap is what remains.** 7.9 µs per pair against TreeDist's 2.6 µs, down from
72.5× to 3.1×. Every pair still re-reads both trees; hoisting that out of the loop, and
running independent pairs in parallel, is the outstanding work. It changes no results.

### Quartet distance

| taxa | quartets | PhyloDistances | Quartet | ratio |
|-----:|---------:|---------------:|--------:|------:|
| 10 | 210 | 9.7 µs | 1.71 ms | 176.4× faster |
| 50 | 230,300 | 821.2 µs | 12.21 ms | 14.9× faster |
| 200 | 64,684,950 | 222.34 ms | 183.46 ms | 1.2× slower |
| 477 | 2,130,031,575 | 7.68 s | 935.14 ms | 8.2× slower |
| 700 | 9,918,641,075 | 37.40 s | — | — |
| 1000 | 41,417,124,750 | 183.30 s | — | — |

The crossover sits near 200 taxa: below it R's per-call overhead dominates and this package
wins comfortably, above it the algorithms decide and the gap widens without bound.

This one is slower than its reference at the sizes that matter, by a margin that grows with
every taxon, and that is expected rather than a defect to be explained away. The
implementation here enumerates all
`binomial(n, 4)` four-taxon subsets and resolves each in constant time from a table of
leaf-to-leaf path lengths — `O(n⁴)`, chosen because it is obviously correct and serves as
the oracle for anything faster. Quartet wraps
[tqDist](https://users-cs.au.dk/cstorm/software/tqdist/), which counts the same quantity
without enumerating it.

So the two columns compare algorithms, not languages, and the gap is the standing argument
for implementing a subquadratic counting scheme. The benchmark exists to say how much that
would be worth.

**Quartet stops at 477 tips.** Above that the quartet count outgrows the signed 32-bit
integers it is accumulated in, and the package refuses with "trees too large for integer
representation". It is worth knowing that the ceiling bites slightly *below* the refusal:
at 477 tips the reported `N` column is already `2 * Q` overflowed to `NA`, while `Q` itself
still fits. Read `Q`, `s`, `d`, `r1`, `r2`, `u`; never `N`. The last rows of the table have
no reference column because of this limit, not because the comparison was skipped.

`run.jl` checks that both implementations returned the same distance for every pair it
timed, and the table says so per row — a benchmark that quietly measured two different
computations would be worthless.
