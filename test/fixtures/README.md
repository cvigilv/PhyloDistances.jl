# Reference values

`rf_quartet.tsv` records the expected [`RobinsonFoulds`](../../src/robinsonfoulds.jl) and
[`QuartetDistance`](../../src/quartet.jl) of a set of named tree pairs. `read.jl` defines
the file's format; `../test_fixtures.jl` reads it and checks every value.

## Columns

| column | meaning |
| --- | --- |
| `case` | a name unique within the file, used to identify a row in a failure message |
| `newick1`, `newick2` | the pair, as Newick strings |
| `rf` | `RobinsonFoulds()(t1, t2)`, an integer |
| `rf_normalized` | `RobinsonFoulds(normalize = true)(t1, t2)`, or `NaN` where the two trees carry no split between them and the divisor is zero |
| `quartet` | `QuartetDistance()(t1, t2)`, an integer |
| `quartet_normalized` | `QuartetDistance(normalize = true)(t1, t2)` |
| `provenance` | the independent derivation of the row's values |

Fields are separated by single tab characters and therefore contain none. Lines beginning
with `#`, and blank lines, are comments; the first line that is neither names the columns.

## Where the values come from

The numeric columns are **generated from the R packages this metric suite reproduces** —
TreeDist for Robinson-Foulds, Quartet for the quartet distance — by
`../../validation/fixture.jl`. That script also re-checks the committed file against them:

```console
$ julia --project=validation validation/fixture.jl           # check, non-zero on any difference
$ julia --project=validation validation/fixture.jl --write    # regenerate
```

Values are compared bitwise rather than to a tolerance, here and in the tests. Both sides
reach them by the same IEEE operations, so a tolerance would hide a real difference rather
than absorb a meaningless one.

The `provenance` column is written by hand and carried through regeneration unchanged. It
states how each row's values follow analytically or from enumerating the quartets, which
is what makes a disagreement between the generator and the file informative rather than
circular: the reference fixes the values, and the derivation says why they are the right
ones. Every row was derived by hand first and agreed with the reference exactly, including
the two six-taxon quartet counts, 8 of 15 and 13 of 15.

The rows exercise, between them, identical trees, a single nearest-neighbor interchange,
trees sharing no split, a rooted input against its unrooted twin, a polytomy against a
resolved tree, branches of length zero, and a pair whose Robinson-Foulds normalizer is
itself zero.

## Portability

Generating the fixture needs R; **reading it does not**, and that is the point of
committing it. The tests that consume this file must pass on any machine with no R, no
network, and no data outside the repository, so the path is resolved relative to the test
file and the trees are inline Newick rather than references to tree files.

The broader randomized comparison against the same two R packages lives in
`../../validation/crosscheck.jl`, which does require R.
