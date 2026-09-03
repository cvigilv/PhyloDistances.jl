# Benchmarks

Compares this package against the R packages whose values it reproduces —
[TreeDist](https://github.com/ms609/TreeDist) for Robinson-Foulds, Jaccard-Robinson-Foulds,
Info-Robinson-Foulds, mutual clustering information and clustering information distance,
and [Quartet](https://github.com/ms609/Quartet) for the quartet distance — so that agreeing
on results does not quietly cost an order of magnitude in speed.

## Running

```console
$ julia --project=benchmark benchmark/run.jl
```

`run.jl` writes Newick files to `benchmark/trees/`, benchmarks this package, invokes
`treedist.R`, `quartet.R`, `jrf.R`, `inforf.R` and `clusteringinfo.R` on the same files,
and renders `results.md`. R is optional — without a working `Rscript` the Julia timings are
reported alone.

**A full run takes well under a minute** on the Julia side; most of the wall clock goes to R,
repeating its calls until the clock's resolution stops mattering.

Neither side's timing includes parsing. `benchmark/trees/` and the `results_*_r.tsv` files
are regenerated on every run and are not tracked; `results.md` is.

### The taxon ramp

```console
$ julia --project=benchmark benchmark/rampbench.jl
```

`rampbench.jl` does the same thing for the two split-matching metrics and the quartet
distance at 16, 64, 256 and 1024 taxa, writing its trees to `benchmark/ramptrees/`, calling
`rampbench.R`, and rendering `ramp.md`. It also leaves both sides' numbers in
`results_ramp_julia.tsv` and `results_ramp_r.tsv`, which are the files to keep a copy of
before a change and compare a later run against; the trees come from a fixed seed, so two
runs measure the same work.

The sizes are powers of two on purpose, which is the whole reason it exists next to
`run.jl` rather than inside it. Anything laid out as an `n × n` table and read at scattered
`[x, y]` — the quartet distance's most-recent-common-ancestor intervals, for one — collides
in the cache when its leading dimension is a power of two, and `run.jl`'s round sizes never
land on one. A ramp that does is what turns that class of problem from invisible into a
ratio that moves.

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
| 10 | 1.6 µs | 46.0 µs | 29.4× faster |
| 50 | 5.7 µs | 53.9 µs | 9.5× faster |
| 200 | 26.4 µs | 109.9 µs | 4.2× faster |
| 1000 | 124.0 µs | 2.20 ms | 17.8× faster |
| all pairs, 40 trees × 60 taxa | 6.62 ms | 1.73 ms | 3.8× slower |

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

**The all-pairs gap is what remains.** 8.5 µs per pair against TreeDist's 2.2 µs, down from
72.5× to 3.8×. Every pair still re-reads both trees; hoisting that out of the loop, and
running independent pairs in parallel, is the outstanding work. It changes no results.

### Quartet distance

| taxa | quartets | PhyloDistances | Quartet | ratio |
|-----:|---------:|---------------:|--------:|------:|
| 10 | 210 | 15.7 µs | 1.38 ms | 88.3× faster |
| 50 | 230,300 | 350.3 µs | 7.23 ms | 20.6× faster |
| 200 | 64,684,950 | 13.16 ms | 102.07 ms | 7.8× faster |
| 477 | 2,130,031,575 | 153.54 ms | 568.43 ms | 3.7× faster |
| 700 | 9,918,641,075 | 449.95 ms | — | — |
| 1000 | 41,417,124,750 | 1.38 s | — | — |
| 1500 | 210,094,780,875 | 4.72 s | — | — |

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

Two costs inside the concordant count were worth more than anything about the algorithm
itself.

**Tabulating a clade's membership is not always worth its `O(n)`.** For each branch point of
the first tree the count builds a prefix sum over the second tree's numbering, which answers
"how many of this clade's leaves lie in that one" in constant time. That pays handsomely at a
branch point that serves thousands of pairs and not at all at a cherry, which scans every
taxon to answer one — and a binary tree has as many cherries as it has deep branch points.
Building the table only where `npairs × |clade|` exceeds `n`, and counting the pairs directly
where it does not, cut about a quarter of the total at every size tested, with no change to
the `O(n³)` bound.

**A power-of-two leading dimension is worth a factor of two to four on its own.** The
most-recent-common-ancestor intervals live in an `n × n` matrix read and written at scattered
`[x, y]`, and when `n` is a power of two those accesses collide in the cache: at 1024 taxa the
count ran four times slower than at 1000, for 2.4% more work. One row of padding removes it.
Packing each interval's two `Int32` endpoints into one `UInt64` — they are always written and
read together — halves the traffic again.

`run.jl` checks that both implementations returned the same distance for every pair it
timed, and the table says so per row — a benchmark that quietly measured two different
computations would be worthless.

### Jaccard-Robinson-Foulds

| taxa | PhyloDistances | TreeDist | ratio |
|-----:|---------------:|---------:|------:|
| 10 | 10.5 µs | 67.0 µs | 6.4× faster |
| 50 | 74.2 µs | 91.9 µs | 1.2× faster |
| 200 | 833.7 µs | 425.9 µs | 2.0× slower |
| 1000 | 22.18 ms | 17.18 ms | 1.3× slower |

An earlier version of this benchmark showed a much worse picture — up to 7× slower at
n=1000, widening rather than narrowing with tree size. Four rounds of work closed most of
that gap, and not by the same route each time; the second is worth being honest about,
because it changed *what* the solver is, not really *how fast* it is.

**Round one fixed real waste.** The score matrix allocated a `BitVector` per cell
(`count(a .& b)` — broadcasting `.&` builds a fresh array before `count` ever runs, once
per pair of splits, `O(n²)` times); replacing it with a direct population count over the
two masks' machine-word chunks (`_countand`) removed the allocation, and hoisting each
split's own marked-taxon count out of the loop removed redundant `O(n)` work that came
with recomputing it per cell. The assignment solver reused its scratch buffers across rows
instead of reallocating them every row. Together these cut allocation at n=1000 from
**245 MB to 10 MB** and roughly halved wall time, for changes that are unambiguous wins —
none of them trade anything away.

**Round two replaced the solver with Jonker & Volgenant's (1987) actual published
algorithm** — the same one TreeDist's C++ implements — rather than the ad hoc row-reduction
heuristic round one had settled on. This was **not primarily a speed change**: on this
package's split-matching workload it lands within a few percent of the heuristic it
replaced, which is itself already close to whatever the assignment solver's inherent cost
is on this kind of matrix. The reason to make the change anyway is that it is the correct,
well-specified algorithm rather than an ad hoc approximation of one — and porting it
surfaced a real bug worth knowing about even if you never touch this code: two rows that
are genuinely tied for the same column's best reduced cost can, after enough potential
updates, present as differing by a single floating-point ulp rather than exactly. Treating
that as "strictly less" sent the pair into an infinite tug-of-war over the column, each
"win" tightening the potential by an amount too small to change anything next time. This
package's own brute-force test suite caught it directly (a hang, not a wrong answer) before
it reached anyone; the fix — a tolerance-scaled comparison (`_jvstrictlyless`) — mirrors a
guard TreeDist's own solver already carries for the same reason, which this project had
initially (wrongly) assumed was specific to TreeDist's integer-quantized costs rather than
a real floating-point degeneracy.

**Round three was not about this metric at all.** The split-hashing change described under
Info-Robinson-Foulds below applies to every metric that builds a split set, and moved this
one from 1.3×/2.6×/2.8× slower to 1.1×/1.9×/2.3× slower at n=50/200/1000 without touching
the solver.

**Round four was about memory layout, in three places, and none of it changed a result.** On
a controlled 16/64/256/1024 ramp it took the metric 1.06×/1.28×/1.44×/1.81× faster and cut
allocation at the top of the ramp from 18.6 MB to 10.8 MB.

- *The score matrix read a thousand separate objects.* Each split is its own `BitVector`, so
  intersecting every pair chased two pointers into unrelated parts of the heap. Copying the
  masks' words into one matrix first — one column per split — makes each pair two contiguous
  runs, and cost microseconds against milliseconds saved.
- *`x^1` is `x`.* The Jaccard score is raised to `k`, and `k = 1` is both the default and the
  only exponent `NyeSimilarity` uses. Choosing the exponent once per matrix rather than
  testing for it once per pair took the matrix from 11.5 ms to 8.5 ms at n=1000.
- *The solver read its cost matrix across the grain.* Every phase of Jonker & Volgenant after
  column reduction scans one row against all columns, which in a Julia `Matrix` strides a
  whole row length per element. An assignment problem and its transpose have the same
  solution with the two sides exchanged, so the solver now reads `cost[j, i]` and swaps its
  two results — the scans run down a column of storage and nothing is copied. A square
  problem, which is what two trees on the same taxa almost always give, also stopped being
  copied into a padded square of the same size.

What is left is the assignment solve itself, which is now most of the metric.

`agree` in `results.md`'s table uses a tolerance rather than exact equality, because the two
implementations solve the underlying assignment problem with independently written
solvers — this package's vendored one over exact `Float64` costs, TreeDist's over costs
quantized for its integer solver — so the two can differ in the last few significant
digits even when both found the true optimum. See `validation/crosscheck.jl`'s
`_closeenough` and its validation report for the correctness case; this section is about
speed only.

### Clustering information

The matching pipeline now removes exact split pairs before it builds the dense score matrix,
and it reads every logarithm in the remaining mutual-information loop from a `log2(0:n)`
table. MCI and CID use the same reduced assignment. CID also uses the table for its tree
entropy sums.

#### MCI before and after

| taxa | before | after | speedup | alloc before | alloc after | TreeDist before / after | ratio before | ratio after |
|-----:|-------:|------:|--------:|-------------:|------------:|------------------------:|-------------:|------------:|
| 10 | 10.8 µs | 9.8 µs | 1.1× | 0.01 MB / 390 | 0.01 MB / 407 | 45.1 / 42.9 µs | 4.2× faster | 4.4× faster |
| 50 | 106.1 µs | 60.9 µs | 1.7× | 0.11 MB / 1,877 | 0.09 MB / 1,897 | 57.0 / 56.0 µs | 1.9× slower | 1.1× slower |
| 200 | 1.33 ms | 271.1 µs | 4.9× | 0.73 MB / 7,309 | 0.41 MB / 7,330 | 214.0 / 211.0 µs | 6.2× slower | 1.3× slower |
| 1000 | 48.66 ms | 2.62 ms | 18.6× | 10.42 MB / 36,161 | 2.89 MB / 36,177 | 5.89 / 5.81 ms | 8.3× slower | 2.2× faster |

#### CID before and after

| taxa | before | after | speedup | alloc before | alloc after | TreeDist before / after | ratio before | ratio after |
|-----:|-------:|------:|--------:|-------------:|------------:|------------------------:|-------------:|------------:|
| 10 | 11.0 µs | 9.5 µs | 1.2× | 0.01 MB / 392 | 0.01 MB / 409 | 101.8 / 98.9 µs | 9.3× faster | 10.4× faster |
| 50 | 105.0 µs | 57.4 µs | 1.8× | 0.11 MB / 1,879 | 0.09 MB / 1,899 | 119.0 / 120.9 µs | 1.1× faster | 2.1× faster |
| 200 | 1.33 ms | 272.6 µs | 4.9× | 0.73 MB / 7,311 | 0.41 MB / 7,332 | 335.9 / 302.1 µs | 4.0× slower | 1.1× faster |
| 1000 | 48.31 ms | 2.64 ms | 18.3× | 10.42 MB / 36,163 | 2.89 MB / 36,179 | 7.64 / 6.19 ms | 6.3× slower | 2.3× faster |

Allocation counts rise by 16 to 21 small objects because the reduced path records unmatched
indices and builds the logarithm table. Allocated bytes fall by 72% at 1000 taxa because the
score matrix and assignment scratch space shrink from 997 by 997 to 193 by 193 on this
seeded pair.

#### What the profile showed

Before editing, the 200-taxon call took 1.36 ms. Split extraction took 83.0 µs, the full
score matrix 746.2 µs, and assignment 342.5 µs. At 1000 taxa the corresponding figures were
48.55 ms end to end, 520.7 µs for split extraction, 24.01 ms for the score matrix, and
21.71 ms for assignment. Taxon indexing added 160.7 µs and 813.2 µs at the two sizes.

The score loop was mostly logarithms. Building only the contingency counts took 98.2 µs at
200 taxa and 8.74 ms at 1000. Scoring precomputed counts took another 576.6 µs and 13.81 ms,
while the same loop with only integer arithmetic took 5.3 µs and 116.0 µs. This ruled out
population counting as the main arithmetic problem.

The two changes were benchmarked separately on pre-extracted splits:

| matching pipeline | 200 taxa | 1000 taxa |
|:------------------|---------:|----------:|
| original | 1.086 ms | 46.83 ms |
| integer-log table only | 597.2 µs | 34.68 ms |
| exact-match removal only | 62.0 µs | 1.72 ms |
| both | 37.0 µs | 1.29 ms |
| both, without redundant bounds checks | 32.6 µs | 1.22 ms |

Exact-match removal is the larger win on these perturbed pairs. They share 152 of 197 splits
at 200 taxa and 804 of 997 at 1000. The logarithm table still matters when few splits match,
and it cuts another quarter from the reduced 1000-taxon pipeline here. The initial profile
also attributed samples to bounds checks in the score loop. Removing checks only where the
matrix axes and contingency-count range prove every access valid cut the reduced matrix by
11% at 1000 taxa and the full call by about 3%.

After these changes, the 1000-taxon call spends about 0.80 ms indexing taxa, 0.52 ms extracting
splits, 0.08 ms identifying exact matches, 0.46 ms constructing the reduced score matrix,
and 0.73 ms solving the reduced assignment. No single local score-loop change dominates any
longer. Further gains require faster shared tree preprocessing or changes to the assignment
solver, both used outside MCI and CID.

Self-comparisons still use the same table-backed entropy expression for the exact score and
the normalizer. Raw MCI therefore equals tree entropy exactly, normalized MCI is exactly
one, and CID is exactly zero. The full 1,140-case cross-check still reports zero mismatches
against TreeDist 2.14.1.

### Info-Robinson-Foulds

| taxa | PhyloDistances | TreeDist | ratio |
|-----:|---------------:|---------:|------:|
| 10 | 10.4 µs | 132.1 µs | 12.7× faster |
| 50 | 58.2 µs | 292.8 µs | 5.0× faster |
| 200 | 257.8 µs | 968.0 µs | 3.8× faster |
| 1000 | 1.43 ms | 6.82 ms | 4.8× faster |

`InfoRobinsonFoulds` weights each split by its phylogenetic information content instead of
counting it as one, but — unlike Jaccard-Robinson-Foulds — matches splits by exact identity,
the same relationship classic Robinson-Foulds already computes as a symmetric-difference
intersection. There is no assignment solve here at all.

It started out 10.0×/2.0× faster at n=10/50 but 1.7×/5.6× *slower* at n=200/1000 — a
steeper falloff than Jaccard-Robinson-Foulds's, despite doing strictly less work per pair.
Two separate `O(n²)` costs were responsible, and neither was the split intersection itself.

**Summing information content recomputed a running sum per split.** `log2rooted` accumulates
`k` logarithms from scratch on every call, which is fine once but not `~n` times: summing a
tree's own split information (`SplitwiseInfo`, computed twice per comparison) cost `O(n²)`
even though nothing else in the metric did. TreeDist tables these values up to 64 tips for
exactly this reason. `SplitInfoTable` computes `log2rooted(0), …, log2rooted(n)` in one pass
and reads each value back in constant time; because it is a prefix scan of the same loop,
adding the same terms in the same order, a tabulated result is *bitwise* identical to the
per-call one rather than merely close. That halved n=1000 and closed the n=200 crossover.

**A split's identity was computed one bit at a time.** What remained was not the arithmetic
at all — after tabling, the information sum was 5 µs of a 20 ms comparison. Splits are
`BitVector`s used as `Dict` keys and `Set` members, and `hash(::BitVector)` walks the vector
element by element: at n=1000 hashing one mask cost **1.383 µs against 16.074 ns for the 16
machine words already backing it**, an 86× difference paid on every lookup. Hashing those
words instead (`SplitKey`) took a comparison from 20.19 ms to 4.77 ms, of which `intersect`
went from 7.67 ms to 78.5 µs — 98× — and building one tree's splits from 5.70 ms to 1.88 ms.

This is sound for the same reason Base's own `==(::BitArray, ::BitArray)` compares `.chunks`
directly: a `BitArray` holds the unused bits of its final word at zero. Length is part of the
identity too, since masks over different taxon counts can hold equal words. The storage is
unchanged — splits are still `BitVector`s — so this carries no n ≤ 64 ceiling, which packing
each split into a single `UInt64` would have.

**Sorting the splits had the same defect the hashing did.** The masks are sorted so that two
trees with the same splits iterate them in the same order, and `isless(::BitVector,
::BitVector)` compares bit by bit: at n=1000 sorting one tree's splits took 1.73 ms of the
1.82 ms the whole `splits` call cost. Ordering on the backing words instead — a different
canonical order, and equally canonical, since the words determine the mask — makes it 0.03 ms.
That is the largest single share of this metric, which it took from 4.77 ms to 1.43 ms at
n=1000, and it applies to every metric that builds a split set.

What remains at n=1000 is the tree walk's per-branch allocation rather than any lookup.

`agree` uses a tolerance rather than exact equality only because TreeDist floors its result
near zero (`.FloorNumericalNoise`), not because either side solves an optimization
differently — see `validation/crosscheck.jl`.
