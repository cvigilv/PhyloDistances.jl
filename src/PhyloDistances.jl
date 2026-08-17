"""
    PhyloDistances

Distance and similarity metrics between phylogenetic trees.

Trees are read from Newick with [NewickTree.jl](https://github.com/arzwa/NewickTree.jl) and
traversed through the [AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl)
interface, so any node type implementing that interface is accepted.

A metric is a value, and [`compare`](@ref) applies it to a pair of trees:

```julia
d = compare(SomeMetric(), tree1, tree2)
```

[`pairwise`](@ref) applies a metric to every pair in a collection.

Metrics are computed under a [`Convention`](@ref), which selects between published
formulations that disagree on normalization or on the treatment of trivial splits.
"""
module PhyloDistances

using AbstractTrees
using NewickTree

export TreeMetric
export Convention, TreeDistConvention, PrimaryConvention
export compare, pairwise

public normalizer, requiresrooted, issimilarity

include("interface.jl")

end
