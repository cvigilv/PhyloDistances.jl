# Agreement with TreeDist

This package reproduces the R package [TreeDist](https://github.com/ms609/TreeDist), so
matching its values is the point and being close is not good enough. `crosscheck.jl`
compares the two on deliberately awkward trees, exactly.

```console
$ julia --project=validation validation/crosscheck.jl [ncases]
```

It writes tree pairs, has `treedist_values.R` compute TreeDist's answers, compares them, and
renders `report.md`. It exits non-zero if any value differs, so it can gate a release.
`julia_cases.tsv` and `r_values.tsv` are regenerated each run and untracked; `report.md` is
kept.

Integers are compared with `==` and floating-point values with `===` — bitwise, no
tolerance. `NaN` must meet `NaN`, which happens when neither tree carries a split and the
normalizer is zero.

## Cases

Random trees alone miss the shapes where implementations diverge, so the fixed part of the
suite covers, at ten sizes from 4 to 400 taxa: star trees, which carry no splits at all;
caterpillars, the most unbalanced shape there is; polytomies on one side and on both, at
three severities; rooted against rooted and rooted against unrooted; and identical trees
against maximally perturbed ones. Random pairs follow.

That mix is not decoration. Counting the clusters of a rooted tree twice — a real bug, since
rooting at a taxon leaves the original root with a single branch below it — passed two
thousand random *unrooted* pairs and was caught by a rooted one.

## Compare in the right direction

R writes its values and Julia parses and compares them. Not the other way around:
`as.numeric` in R does not reliably round-trip R's own `%.17g` output, landing half an ulp
away on values such as `92/94`. Feeding Julia's output through it reports differences that
do not exist — and would equally hide real ones. Julia's parser round-trips correctly, so
the comparison happens there.
