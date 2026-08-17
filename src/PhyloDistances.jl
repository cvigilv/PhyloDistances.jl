"""
    PhyloDistances

Distance and similarity metrics between phylogenetic trees.

Trees are read from Newick with [NewickTree.jl](https://github.com/arzwa/NewickTree.jl) and
traversed through the [AbstractTrees.jl](https://github.com/JuliaCollections/AbstractTrees.jl)
interface, so any node type implementing that interface is accepted.
"""
module PhyloDistances

using AbstractTrees
using NewickTree

end
