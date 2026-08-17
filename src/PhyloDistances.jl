"""
    PhyloDistances

Distance and similarity metrics between phylogenetic trees.

Trees are read from Newick with [NewickTree.jl](https://github.com/arzwa/NewickTree.jl) and
traversed through the [AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl)
interface, so any node type implementing that interface is accepted.

Metrics implement the [Distances.jl](https://github.com/JuliaStats/Distances.jl) interface.
A metric is a value carrying everything that selects a variant of the computation, and
applying it to two trees compares them:

```julia
RobinsonFoulds()(tree1, tree2)
pairwise(RobinsonFoulds(), trees)
```

`pairwise`, `colwise`, `evaluate` and `result_type` are re-exported from Distances.jl, so
loading both packages is safe — the names refer to the same functions.

Distances subtype [`TreeMetric`](@ref); quantities that are largest on identical trees
subtype [`TreeSimilarity`](@ref) instead. Both are computed under a [`Convention`](@ref),
which selects between published formulations that disagree on normalization or on the
treatment of trivial splits.
"""
module PhyloDistances

using AbstractTrees
using Distances: Distances, colwise, colwise!, evaluate, pairwise, pairwise!, result_type
using NewickTree

export TreeMetric, TreeSimilarity, TreeComparison
export Convention, TreeDistConvention, PrimaryConvention

# Distances.jl is the interface these metrics implement; re-exporting means a user loading
# both packages sees one set of functions rather than a name clash.
export colwise, colwise!, evaluate, pairwise, pairwise!, result_type

public convention, isnormalized, normalizer, requiresrooted, issimilarity

include("interface.jl")

end
