# Benchmarks

Compares this package against [TreeDist](https://github.com/ms609/TreeDist), the R package
whose values it reproduces, so that agreeing on results does not quietly cost an order of
magnitude in speed.

## Running

```console
$ julia --project=benchmark benchmark/run.jl
```

`run.jl` writes Newick files to `benchmark/trees/`, benchmarks this package, invokes
`treedist.R` on the same files, and renders `results.md`. TreeDist is optional — without a
working `Rscript` the Julia timings are reported alone.

Neither side's timing includes parsing. `benchmark/trees/` and `results_r.tsv` are
regenerated on every run and are not tracked; `results.md` is.

Memory is measured differently on each side and the two figures are **not** comparable.
Julia reports bytes allocated by the call; R reports the peak its garbage collector saw,
which includes everything already resident.

## Reading the results

The picture at the time of writing, from `results.md`:

| taxa | PhyloDistances | TreeDist | ratio |
|-----:|---------------:|---------:|------:|
| 10 | 1.5 µs | 37.2 µs | 25× faster |
| 50 | 5.5 µs | 46.0 µs | 8× faster |
| 200 | 26.2 µs | 115.9 µs | 4× faster |
| 1000 | 120 µs | 2.41 ms | 20× faster |
| all pairs, 40 trees × 60 taxa | 6.27 ms | 1.73 ms | 3.6× slower |

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

**The all-pairs gap is what remains.** 8.0 µs per pair against TreeDist's 2.2 µs, down from
72.5× to 3.6×. Every pair still re-reads both trees; hoisting that out of the loop, and
running independent pairs in parallel, is the outstanding work. It changes no results.

None of this affects correctness: values agree with TreeDist across hand-checked pairs and a
randomized comparison over 1,500 trees, half of them rooted, and separately against the
split-set formulation the encoding replaced.
