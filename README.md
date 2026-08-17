# PhyloDistances.jl

Distance and similarity metrics between phylogenetic trees.

Trees are read from Newick with [NewickTree.jl](https://github.com/arzwa/NewickTree.jl) and
traversed through the
[AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl) interface, so any
node type implementing that interface is accepted.

## Interface

A metric is a value, and `compare` applies it to two trees:

```julia
using PhyloDistances

d = compare(SomeMetric(), tree1, tree2)
```

Parameters that select a variant of a metric are fields of the metric type, so a
constructed metric specifies the computation completely:

```julia
compare(SomeParameterizedMetric(k = 2), tree1, tree2)
```

`pairwise` applies a metric to every pair in a collection, returning a square matrix that
shares the axes of its input:

```julia
D = pairwise(SomeMetric(), trees)
```

### Conventions

Several tree metrics have more than one formulation in the literature, differing in
normalization, in the treatment of trivial splits, or in the base of the logarithm. The
`convention` keyword selects which is computed, so a value can be traced back to a specific
definition:

```julia
compare(SomeMetric(), tree1, tree2; convention = :treedist)  # default
compare(SomeMetric(), tree1, tree2; convention = :primary)
```

`:treedist` reproduces the R package [TreeDist](https://github.com/ms609/TreeDist) and is
the default; `:primary` follows the source that first defined the metric, where the two
differ.

### Normalization

`normalize = true` divides the result by the largest value the metric can take on trees
with the given taxon set, placing it on a scale comparable across taxon set sizes:

```julia
compare(SomeMetric(), tree1, tree2; normalize = true)
```

### Traits

`PhyloDistances.issimilarity(metric)` reports whether larger values mean *more* similar
trees; `PhyloDistances.requiresrooted(metric)` reports whether the metric is defined on
rooted trees.
