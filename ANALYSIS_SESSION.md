# Session Handoff — 2026-08-19

## Project maturity target

`package` — PhyloDistances.jl (new package, does not extend an existing one)

## What was just completed

CHUNK-008: validate-rf-and-quartet — the "usable today" checkpoint, which closes the
critical path. Robinson-Foulds and quartet distance are now callable and validated.

`test/fixtures/rf_quartet.tsv` is a tab-separated table of nine named tree pairs with the
expected `RobinsonFoulds` and `QuartetDistance` value of each, raw and normalized, plus a
per-row provenance note. Its numeric columns are **generated from TreeDist 2.14.1 and
Quartet 1.3.0** by the new `validation/fixture.jl`, which also re-checks the committed file
against them and exits non-zero on any difference. `test/fixtures/read.jl` holds the one
reader both sides use; `test/test_fixtures.jl` checks every value and then the all-pairs API
over a five-tree collection.

## Key decisions made

- **The reference fixes the values; the hand derivation says why they are the right ones.**
  Every row was derived analytically or by enumerating its quartets before anything was run,
  and the file is then generated from TreeDist and Quartet. The two agreed bitwise on all
  nine rows. Either alone would be weaker: values read off the implementation would preserve
  a bug happily, and values taken from the reference alone would not catch reading the
  reference wrongly — which counts `QuartetStatus` reports to combine is itself a choice.
- **Only the numeric columns are generated.** `provenance` is written by hand and carried
  through regeneration unchanged, so `--write` can never quietly rewrite the justification
  along with the number.
- **Plain TSV, parsed by `split(line, '\t')`**, rather than a serialization format. It adds
  no dependency, diffs cleanly, and is readable in the file. Lines starting with `#` are
  comments and the first other line names the columns, which the parser asserts.
- **Everything is compared bitwise, in the tests as in `validation/`.** Both sides reach
  these values by the same IEEE operations, so a tolerance would hide a real difference
  rather than absorb a meaningless one. `NaN` agrees only with `NaN`.
- **`NaN` is a fixture value, not a special case.** The two-stars row has no split in either
  tree, so its RF normalizer is zero; the parser reads `NaN` and the comparison helper
  routes it to `isnan`.

## The nine cases

| case | covers |
| --- | --- |
| `identical-binary` | a tree against itself |
| `four-taxa-one-quartet` | the smallest pair that can differ |
| `one-nni-apart` | a single nearest-neighbor interchange |
| `star-vs-resolved` | maximal quartet distance, `C(n,4)` |
| `rooted-vs-unrooted-twin` | the rooted-input warning path |
| `polytomy-vs-resolved` | a polytomy against a resolved tree |
| `no-shared-splits` | two caterpillars with nothing in common |
| `zero-length-branches` | lengths ignored by both metrics |
| `both-stars` | a zero RF normalizer, giving `NaN` |

## Hand-deriving a quartet count

Worth reusing, since it is what let the fixture be independent. On an unrooted
**caterpillar**, number each taxon by the internal node it hangs from along the path; a
quartet's topology then pairs the two lowest-numbered taxa against the two highest. So
`(A,B,(((C,D),E),F));` numbers A→1, B→1, F→2, E→3, C→4, D→4, and the quartet `{B,D,E,F}` is
read straight off as `BF|ED`. All fifteen quartets of a six-taxon pair take a couple of
minutes this way, against arithmetic on three path-length sums per quartet.

Where a tree has a polytomy the caterpillar trick does not apply and the split-set rule
does: a tree resolves a quartet exactly when one of its splits cuts it two against two.

## State of the codebase

- Files created: `test/fixtures/rf_quartet.tsv`, `test/fixtures/read.jl`,
  `test/fixtures/README.md`, `test/test_fixtures.jl`, `validation/fixture.jl`
- Files modified: `test/runtests.jl`, `validation/README.md`, `validation/report.md`
- Files deleted: `test/fixtures/.gitkeep` (the directory now has real content)
- Package loads cleanly: yes, on Julia 1.12.6
- Test suite passes: yes — **1274 tests**, up from 963, with no R and no network
- R crosscheck passes: yes — `validation/crosscheck.jl` over **1140 cases, zero mismatches**,
  and `validation/fixture.jl` confirms all nine fixture rows
- Entry points: `RobinsonFoulds()(t1, t2)`, `QuartetDistance()(t1, t2)`,
  `pairwise(metric, trees)`
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

## Watch out for

- **The R references need the Nix R 4.4.2, not the R on `PATH`.** `/usr/local/bin/R` is a
  CRAN R 4.6 whose library has neither TreeDist nor Quartet. Both are installed for R 4.4
  and load only under the Nix build, so prefix every validation command:

  ```console
  $ PATH=/nix/store/jaqvbj23b52yl0qgcrrb4ysbxdlqlbv5-R-4.4.2-wrapper/bin:$PATH \
      julia --project=validation validation/crosscheck.jl
  ```

  Two ways this misleads. An R that cannot see a package reports it as **not installed**
  rather than as installed elsewhere, so `packageVersion` failing proves nothing about the
  machine — check `.libPaths()` and `R_LIBS_USER`, which is version-specific. And loading a
  package built for a different R **aborts the session** ("An irrecoverable exception
  occurred") instead of raising a catchable error.
- **`readnw` does not accept a `SubString`.** `split(line, '\t')` yields `SubString{String}`
  and `readnw` fails on one with a misleading "no trailing semicolon?" message wrapping a
  `MethodError`. `test/fixtures/read.jl` converts with `String(...)`; any other code reading
  Newick out of a parsed file must do the same.
- **The test suite prints ~400 uncaptured rooting warnings** from randomized cases in
  `test/test_robinsonfoulds.jl`. Correct behavior, unreadable output, and a genuine warning
  would be lost in it. Recorded in Open Questions for CHUNK-023.
- **`splits` excludes trivial splits by default**; pass `trivial = true` to sum over all
  branches.
- **A node with one branch below it inflates path lengths across it** — a rooted tree's root
  is one. Harmless for quartets, but it silently inflated every RF distance involving a
  rooted input in CHUNK-006. Treat each new traversal on its own terms.
- **Randomized checks must include rooted and polytomous trees.** `randomtree` produces
  neither; 2,000 random unrooted pairs missed the CHUNK-006 root bug that seven hand-picked
  pairs caught.
- **Derive normalizers from the reference, never from the literature** — TreeDist normalizes
  RF by `n1 + n2`, not `2(n-3)`.
- **`Bool <: Real` in Julia**, so `normalize = true` means the metric's own scheme and
  `normalize = 1` means divide by one.
- **Length arithmetic across rootings is inexact**, so branch-length metrics assert with `≈`.
- **`NewickTree.Node(id, data)` leaves `parent` and `children` undefined**; build with
  `push!(parent, child)`.
- **Do not reintroduce a bare `using AbstractTrees` or `using NewickTree`** — both export
  `children`, `getroot`, `isroot` and `print_tree`.
- **Julia 1.12 is the declared floor; there is no LTS back-support.**

## Next chunk

CHUNK-009: branch-score-distance — the Kuhner-Felsenstein branch score, the Euclidean
distance between the two trees' branch-length vectors indexed by split, with an absent split
contributing length zero.

`WeightedRobinsonFoulds` in `src/robinsonfoulds.jl` is the L1 form of exactly this and is
already written; the branch score is its L2 counterpart, so the chunk is largely a second
reduction over `_weightedsplits`. Note what that existing code establishes: the sum runs
over the **union** of the two split sets, `trivial = true` so pendant branches count, and a
tree without branch lengths is rejected rather than yielding `NaN`.

**Neither metric has a reference implementation available.** TreeDist has no
branch-length-weighted distances at all — they belong to phangorn (`wRF.dist`, `KF.dist`).
phangorn is **not** installed: checked under the Nix R 4.4.2 against the library that does
hold TreeDist 2.14.1, Quartet 1.3.0 and TreeTools 2.4.0. Both conventions therefore
coincide, and validation is by hand: identical lengths give zero, one edge differing by δ
gives δ, and the metric axioms hold on random triples.

Installing phangorn there would give the branch-length family a reference, and the recipe
that installed TreeDist and Quartet is known to work — that is CHUNK-032, and doing it
before CHUNK-009 rather than after would let the branch score be validated the same way
Robinson-Foulds and the quartet distance now are, through `validation/`.
