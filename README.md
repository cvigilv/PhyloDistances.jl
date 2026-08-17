# PhyloDistances.jl

Distance and similarity metrics between phylogenetic trees.

Trees are read from Newick with [NewickTree.jl](https://github.com/arzwa/NewickTree.jl) and
traversed through the
[AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl) interface, so any
node type implementing that interface is accepted.

## Interface

Metrics implement the [Distances.jl](https://github.com/JuliaStats/Distances.jl) interface.
A metric is a value, and applying it to two trees compares them:

```julia
using PhyloDistances

d = SomeMetric()(tree1, tree2)
d = evaluate(SomeMetric(), tree1, tree2)   # the same thing
```

Everything that selects a variant of the computation is a field of the metric, so a
constructed metric specifies the computation completely and no keyword arguments are needed
when applying it:

```julia
SomeMetric(; convention = :primary, normalize = true)(tree1, tree2)
SomeParameterizedMetric(2; normalize = true)(tree1, tree2)
```

That is what lets the whole Distances.jl toolkit work unchanged:

```julia
D = pairwise(SomeMetric(), trees)          # all-pairs matrix
pairwise!(D, SomeMetric(), trees)          # in place
colwise(SomeMetric(), trees, references)   # elementwise along two collections
```

`pairwise`, `pairwise!`, `colwise`, `colwise!`, `evaluate` and `result_type` are
re-exported from Distances.jl, so loading both packages is safe: the names refer to the
same functions rather than clashing.

Distance matrices produced this way feed directly into ecosystem tools that expect one —
clustering, multidimensional scaling, and so on.

### Distances and similarities

Distances subtype `TreeMetric`, which is a `Distances.SemiMetric`.

Quantities that are *largest* on identical trees — mutual clustering information, shared
phylogenetic information, Nye similarity, maximum agreement subtree size — subtype
`TreeSimilarity` instead, which deliberately sits outside the Distances.jl hierarchy.
`Distances.PreMetric` requires `d(x, x) == 0`, which a similarity does not satisfy, and
declaring one a `SemiMetric` would be actively wrong: `pairwise` takes the zero diagonal on
faith rather than computing it, so the result would report zero self-similarity. Both kinds
are used identically; only the type hierarchy differs.

### Conventions

Several tree metrics have more than one formulation in the literature, differing in
normalization, in the treatment of trivial splits, or in the base of the logarithm. The
`convention` field selects which is computed, so a value can be traced back to a specific
definition:

```julia
SomeMetric(; convention = :treedist)   # default
SomeMetric(; convention = :primary)
```

`:treedist` reproduces the R package [TreeDist](https://github.com/ms609/TreeDist);
`:primary` follows the source that first defined the metric, where the two differ.

### Normalization

`normalize = true` divides the result by the largest value the metric can take on trees
with the given taxon set, placing it on a scale comparable across taxon set sizes:

```julia
SomeMetric(; normalize = true)(tree1, tree2)
```

### Traits

`PhyloDistances.issimilarity(m)` reports whether larger values mean *more* similar trees;
`PhyloDistances.requiresrooted(m)` reports whether the metric is defined on rooted trees;
`PhyloDistances.convention(m)` and `PhyloDistances.isnormalized(m)` report how a
constructed metric is configured.
