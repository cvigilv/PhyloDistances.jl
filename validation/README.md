# Agreement with TreeDist and Quartet

This package reproduces two R packages — [TreeDist](https://github.com/ms609/TreeDist) for
Robinson-Foulds and [Quartet](https://github.com/ms609/Quartet) for the quartet distance —
so matching their values is the point and being close is not good enough. `crosscheck.jl`
compares them on deliberately awkward trees, exactly.

```console
$ julia --project=validation validation/crosscheck.jl [ncases]
```

It generates tree pairs, computes every quantity on both sides, and renders `report.md`. It
exits non-zero if any value differs, so it can gate a release.

## The committed fixture is generated from the same references

`../test/fixtures/rf_quartet.tsv` holds the expected values for a small set of named tree
pairs, and the portable test suite reads it with no R installed. Those values are not
invented: `fixture.jl` computes each one with TreeDist and Quartet and either verifies the
committed file against them or rewrites it.

```console
$ julia --project=validation validation/fixture.jl           # check, non-zero on any difference
$ julia --project=validation validation/fixture.jl --write    # regenerate
```

The two scripts divide the work between them. `crosscheck.jl` sweeps hundreds of generated
pairs and reports; `fixture.jl` pins a handful of hand-chosen ones into a file that travels
with the tests. Only the numeric columns are generated — the fixture's `provenance` column
is written by hand and says how each row follows analytically, so a disagreement between
generator and file means one of them is wrong rather than that a value has merely moved.

## R must be the one the packages were built for

TreeDist and Quartet load compiled code, which R will not accept across a version change:
a package built for R 4.4 aborts the session when loaded under 4.6. Both scripts take R
from `PATH` and set `R_HOME` from it, so putting the right R first is enough:

```console
$ PATH=/path/to/the/right/R/bin:$PATH julia --project=validation validation/fixture.jl
```

Check with `Rscript -e 'print(.libPaths())'` that the library holding TreeDist is one the
running R will search — `R_LIBS_USER` is version-specific, and an R that cannot see the
package reports it as not installed rather than as built for another version.

Integers are compared with `==` and floating-point values with `===` — bitwise, no
tolerance. `NaN` must meet `NaN`, which happens when neither tree carries a split and the
normalizer is zero.

## Values cross in memory, not as text

R is called through [RCall](https://github.com/JuliaInterop/RCall.jl), so every value
arrives as a machine double or integer and no formatting is involved. This is not merely
tidier. R's `as.numeric` does not reliably round-trip R's own `%.17g` output — it can land
half an ulp away on values such as `92/94` — so a comparison routed through a file has to be
written in one particular direction (R writes, Julia parses) to avoid reporting differences
that do not exist, and to avoid hiding real ones. Passing values in memory removes the
question rather than working around it.

Two further conveniences follow: R's `stop()` surfaces as a catchable Julia exception rather
than a subprocess exit code, and `NA_integer_` arrives as `missing` rather than the string
`"NA"`.

RCall records where R lives when it is built, which on NixOS goes stale whenever a rebuild
moves the store path. `crosscheck.jl` sets `R_HOME` from the `R` on `PATH` before loading
RCall, so the two stay in step without a rebuild.

**The benchmarks deliberately do not use RCall** — see `benchmark/README.md`. An RCall round
trip costs 12–30 µs against roughly 1 µs for the same call timed inside R, which is
comparable to the fastest quantities being benchmarked. Timing has to happen inside R;
comparing values does not.

## The two R packages collide

TreeDist and Quartet both export a function named `RobinsonFoulds`, meaning different
things, and whichever is attached second wins. Quartet's is deprecated and expects a quartet
status vector, so calling it with two trees fails somewhere deep inside with
`subscript out of bounds`. Every call in `crosscheck.jl` is namespace-qualified
(`TreeDist::RobinsonFoulds`, `Quartet::QuartetStatus`) so that load order cannot decide
which function runs.

## What the quartet reference actually reports

`Quartet::QuartetStatus` does not return a distance. It returns counts over the four-taxon
subsets: `s` resolved the same way by both trees, `d` resolved in conflicting ways, `r1` and
`r2` resolved by only the first or only the second tree, `u` unresolved in both, out of `Q`
in total.

`QuartetDistance` counts a subset that only one tree resolves as a difference, so the
reference value is **`d + r1 + r2`**. The crosscheck also asserts `Q == binomial(n, 4)`, so
the two sides are held to agree on the divisor as well as on the distance.

The `N` column is never read. It is `2 * Q`, which overflows R's 32-bit integers at 477 tips
while `Q` itself still fits.

## Cases

Random trees alone miss the shapes where implementations diverge, so the fixed part of the
suite covers, at ten sizes from 4 to 400 taxa: star trees, which carry no splits at all;
caterpillars, the most unbalanced shape there is; polytomies on one side and on both, at
three severities; rooted against rooted and rooted against unrooted; and identical trees
against maximally perturbed ones. Random pairs follow.

That mix is not decoration. Counting the clusters of a rooted tree twice — a real bug, since
rooting at a taxon leaves the original root with a single branch below it — passed two
thousand random *unrooted* pairs and was caught by a rooted one.

Polytomies matter for the quartet distance in particular, because that is the only place the
two published counting conventions disagree: on binary trees nothing is unresolved, so
`r1 = r2 = u = 0` and every definition coincides.
