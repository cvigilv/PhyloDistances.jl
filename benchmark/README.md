# Benchmarks

Compares this package against the R packages whose values it reproduces —
[TreeDist](https://github.com/ms609/TreeDist) for Robinson-Foulds and
Jaccard-Robinson-Foulds, and [Quartet](https://github.com/ms609/Quartet) for the quartet
distance — so that agreeing on results does not quietly cost an order of magnitude in
speed.

## Running

```console
$ julia --project=benchmark benchmark/run.jl
```

`run.jl` writes Newick files to `benchmark/trees/`, benchmarks this package, invokes
`treedist.R` and `quartet.R` on the same files, and renders `results.md`. R is optional —
without a working `Rscript` the Julia timings are reported alone.

**A full run takes well under a minute** on the Julia side; most of the wall clock goes to R,
repeating its calls until the clock's resolution stops mattering.

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
| 10 | 210 | 17.6 µs | 1.40 ms | 79.7× faster |
| 50 | 230,300 | 427.2 µs | 11.25 ms | 26.3× faster |
| 200 | 64,684,950 | 17.24 ms | 160.91 ms | 9.3× faster |
| 477 | 2,130,031,575 | 218.29 ms | 899.19 ms | 4.1× faster |
| 700 | 9,918,641,075 | 681.28 ms | — | — |
| 1000 | 41,417,124,750 | 1.96 s | — | — |
| 1500 | 210,094,780,875 | 7.33 s | — | — |

`QuartetDistance` defaults to `algorithm = :fast`, an `O(n³)` scheme that counts concordant
quartets by reducing each to a rooted triple under every possible outgroup, rather than
enumerating all `binomial(n, 4)` subsets directly (`algorithm = :naive`, still `O(n⁴)`,
remains available as the correctness oracle the fast path is tested against). Quartet wraps
[tqDist](https://users-cs.au.dk/cstorm/software/tqdist/), which counts the same quantity in
`O(n log n)` without enumerating it either — the faster of two subquadratic algorithms, not
a subquadratic one against a quartic one, and the remaining gap past 477 taxa reflects that
rather than a defect here. Where the two overlap, `:fast` now wins outright rather than only
below R's per-call-overhead crossover, as the `O(n⁴)` enumeration used to.

**Quartet stops at 477 tips.** Above that the quartet count outgrows the signed 32-bit
integers it is accumulated in, and the package refuses with "trees too large for integer
representation". It is worth knowing that the ceiling bites slightly *below* the refusal:
at 477 tips the reported `N` column is already `2 * Q` overflowed to `NA`, while `Q` itself
still fits. Read `Q`, `s`, `d`, `r1`, `r2`, `u`; never `N`. The rows past that limit have
no reference column because of it, not because the comparison was skipped.

`run.jl` checks that both implementations returned the same distance for every pair it
timed, and the table says so per row — a benchmark that quietly measured two different
computations would be worthless.

### Jaccard-Robinson-Foulds

| taxa | PhyloDistances | TreeDist | ratio |
|-----:|---------------:|---------:|------:|
| 10 | 11.0 µs | 56.5 µs | 5.1× faster |
| 50 | 110.5 µs | 82.0 µs | 1.3× slower |
| 200 | 1.29 ms | 465.2 µs | 2.8× slower |
| 1000 | 56.1 ms | 17.7 ms | 3.2× slower |

An earlier version of this benchmark showed a much worse picture — up to 7× slower at
n=1000, widening rather than narrowing with tree size. Two fixes closed most of that gap.

**The score matrix was allocating a `BitVector` per cell.** `JaccardRobinsonFoulds` scores
every pair of the two trees' splits, an `O(n²)` matrix in the split count; building it
called `count(a .& b)` per cell, and broadcasting `.&` allocates a fresh `BitVector` to hold
the result before `count` ever runs. Replacing that with a direct population count over the
two masks' machine-word chunks (`_countand` in `src/generalizedrf.jl`) removed the
allocation entirely, and hoisting each split's own marked-taxon count out of the `O(n²)`
loop — it depends on the row or column alone, never on the pair — removed the redundant
`O(n)` work that came with recomputing it per cell. Together with reusing the assignment
solver's scratch buffers across rows instead of reallocating them (`src/assignment.jl`),
this cut allocation at n=1000 from **245 MB to 10 MB**.

**The assignment solver itself was doing far more search than necessary.** Its potentials
started at zero, so the first attempt at a match for any row began from no information at
all. Seeding the row potentials from each row's own minimum cost, and processing rows in
ascending order of that minimum — both valid regardless of starting point or order, since
neither changes which matching is optimal, only how much work finding it takes — cut the
number of augmenting-path steps taken on a ~1000-split matching by roughly half. (Column
reduction, the natural next step for the *square* assignment problem, turned out to be
unsound here once rows and columns are unequal in count — see the comment above the row
reduction in `src/assignment.jl` for the case that caught it.)

What is left is close to the assignment solver's inherent cost on this problem: TreeDist's
C++ solver is genuinely faster per comparison it makes, and closing the remaining gap
further would mean porting more of Jonker & Volgenant's (1987) preprocessing — the
"augmenting row reduction" step, not attempted here — rather than removing more overhead.

`agree` in `results.md`'s table uses a tolerance rather than exact equality, because the two
implementations solve the underlying assignment problem with independently written
solvers — this package's vendored one over exact `Float64` costs, TreeDist's over costs
quantized for its integer solver — so the two can differ in the last few significant
digits even when both found the true optimum. See `validation/crosscheck.jl`'s
`_closeenough` and its validation report for the correctness case; this section is about
speed only.
