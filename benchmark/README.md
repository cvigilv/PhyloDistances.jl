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

| taxa | ratio to TreeDist |
|-----:|------------------:|
| 10 | 0.4× (faster) |
| 50 | 2.5× |
| 200 | 10.8× |
| 1000 | 9.9× |
| all pairs, 40 trees × 60 taxa | 72.5× |

Three things account for the shape of this.

**Small trees favour Julia.** At ten taxa the work is trivial and R's per-call overhead
dominates, so the comparison measures interpreter startup rather than algorithms.

**Set operations dominate a single comparison.** Splitting the work at 200 taxa gives
roughly 180 µs to build the taxon index, 500 µs to extract both split sets, and 600 µs in
`symdiff` — nearly half the total. Splits are `BitVector`s used as hash keys, so each
hash and comparison costs O(n) and the set operations come to O(n²) overall. Packing a
split into a `UInt64` when there are at most 64 taxa would remove that factor entirely, and
the `Splits` API hides the representation, so it can change without touching callers.

**The all-pairs gap is a different problem.** 179 µs per pair against TreeDist's 2.5 µs is
far worse than the single-pair ratio, because every pair rebuilds both trees' split sets
from scratch: m² extractions where m would do. Hoisting that preprocessing out of the loop
is the single largest available win, and it changes no results — only how often the same
work is repeated. TreeDist additionally uses the Day (1985) cluster-table algorithm with a
batched C++ path, which is a deeper difference.

None of this affects correctness: Robinson-Foulds agrees with TreeDist across hand-checked
pairs and a randomized comparison over several hundred trees.
