"""
    PhyloDistances

Distance and similarity metrics between phylogenetic trees.

Trees are read from Newick with [NewickTree.jl](https://github.com/arzwa/NewickTree.jl) and
traversed through the [AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl)
interface, so any node type implementing that interface is accepted once
[`PhyloDistances.taxonlabel`](@ref) is defined for it.

Metrics implement the [Distances.jl](https://github.com/JuliaStats/Distances.jl) interface.
A metric is a value carrying everything that selects a variant of the computation, and
applying it to two trees compares them:

```julia
RobinsonFoulds()(tree1, tree2)
pairwise(RobinsonFoulds(), trees)
```

`pairwise`, `colwise`, `evaluate` and `result_type` are re-exported from Distances.jl, and
`readnw` from NewickTree.jl, so loading those packages alongside this one is safe — the
names refer to the same functions.

Distances subtype [`TreeMetric`](@ref); quantities that are largest on identical trees
subtype [`TreeSimilarity`](@ref) instead. Both are computed under a [`Convention`](@ref),
which selects between published formulations that disagree on normalization or on the
treatment of trivial splits.

Trees are unrooted unless a metric declares otherwise via
[`PhyloDistances.requiresrooted`](@ref); see [`PhyloDistances.isrooted`](@ref) for how
rooting is read from a tree and what happens when it does not match.
"""
module PhyloDistances

# AbstractTrees and NewickTree both export `children`, `getroot`, `isroot` and
# `print_tree`; importing the modules rather than their names keeps every use unambiguous.
using AbstractTrees: AbstractTrees
using Distances: Distances, colwise, colwise!, evaluate, pairwise, pairwise!, result_type
using NewickTree: NewickTree, readnw
using Random: Random

export TreeMetric, TreeSimilarity, TreeComparison
export Convention, TreeDistConvention, PrimaryConvention
export TaxonIndex, taxa, taxonindex
export Splits, splits
export randomtree, perturb
export RobinsonFoulds, WeightedRobinsonFoulds, InfoRobinsonFoulds
export QuartetDistance
export NyeSimilarity, JaccardRobinsonFoulds
export MutualClusteringInfo, ClusteringInfoDistance

# Distances.jl is the interface these metrics implement, and NewickTree.jl reads the trees
# they consume; re-exporting means a user loading those packages too sees one set of
# functions rather than a name clash.
export colwise, colwise!, evaluate, pairwise, pairwise!, result_type
export readnw

public branchlength, clusteringentropy, convention, incidencematrix, isnormalized,
    isrooted, issimilarity, istrivial, jointentropy, log2rooted, log2unrooted,
    mutualinformation, nni!, normalization, normalizer, normalizerinfo, requiresrooted,
    SplitInfoTable, splitinfo, splitmatching, taxonlabel, taxonlabels

include("taxa.jl")
include("random.jl")
include("splits.jl")
include("information.jl")
include("clustertable.jl")
include("assignment.jl")
include("interface.jl")
include("splitmatching.jl")
include("robinsonfoulds.jl")
include("quartet.jl")
include("generalizedrf.jl")
include("clusteringinformation.jl")

end
