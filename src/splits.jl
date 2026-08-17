"""
    branchlength(node) -> Float64

The length of the branch immediately above `node`, or `NaN` where the tree records none.

Branch lengths are not part of the AbstractTrees.jl interface, so this is the hook a node
type implements alongside [`taxonlabel`](@ref) to be usable here.
"""
function branchlength(node)
    throw(ArgumentError(
        "PhyloDistances.branchlength is not defined for $(typeof(node)); " *
        "define it to return the length of the branch above a node"
    ))
end

branchlength(node::NewickTree.Node) = NewickTree.distance(node)

"""
    Splits{L}

The bipartitions a tree induces on its taxa, with the branch length supporting each.

Every branch splits the taxa in two: those below it and the rest. A split is stored as a
`BitVector` over the positions of a [`TaxonIndex`](@ref), oriented canonically so that the
first taxon is never a member. Two trees' splits are therefore comparable exactly when
they were built against the same index, which the set operations enforce.

Orientation makes the representation inherently **unrooted**: a rooted tree's two
root branches describe the same bipartition and collapse to a single split whose length is
their sum, reconstructing the branch an unrooted reading would see. Rooted metrics take
their information from path lengths rather than from splits.

Indexing with a mask gives the supporting branch length; iterating gives the masks.
Construct with [`splits`](@ref).
"""
struct Splits{L}
    index::TaxonIndex{L}
    masks::Vector{BitVector}
    lengths::Dict{BitVector,Float64}

    Splits{L}(index, masks, lengths) where {L} = new{L}(index, masks, lengths)
end

"""
    istrivial(mask, ntaxa) -> Bool

Whether a canonically oriented `mask` separates a single taxon from the rest.

Trivial splits correspond to the pendant branches leading to each leaf. Every tree on a
given taxon set has all of them, so they distinguish no two trees and carry no topological
information — though their lengths still matter to metrics built on branch lengths.
"""
function istrivial(mask::AbstractVector{Bool}, ntaxa::Integer)
    marked = count(mask)
    return marked == 1 || marked == ntaxa - 1
end

# Orienting every split away from the first taxon gives one representative per
# bipartition, so the two branches either side of a root reduce to the same split.
_canonical(mask::BitVector) = first(mask) ? .!mask : copy(mask)

# Walks the tree bottom-up, handing `record!` the split and length of the branch above
# every node except the root, which has no branch above it. Returns the taxa below `node`.
function _walksplits!(record!, node, index::TaxonIndex, ntaxa::Integer)
    kids = AbstractTrees.children(node)

    if isempty(kids)
        mask = falses(ntaxa)
        mask[index[taxonlabel(node)]] = true
        return mask
    end

    mask = falses(ntaxa)
    for child in kids
        below = _walksplits!(record!, child, index, ntaxa)
        mask .|= below
        record!(below, branchlength(child))
    end
    return mask
end

"""
    splits(tree; trivial = false) -> Splits
    splits(tree, index::TaxonIndex; trivial = false) -> Splits

The bipartitions `tree` induces on its taxa.

Supplying an `index` builds the splits against that taxon ordering, which is how two trees
are made comparable; omitting it derives the ordering from `tree` alone.

Trivial splits — the pendant branches separating one taxon from the rest — are excluded by
default, since every tree on the same taxa has all of them. Pass `trivial = true` to keep
them, as metrics that compare branch lengths rather than topologies need.

Branches describing the same bipartition contribute one split whose length is their sum.
This is what a rooted tree's two root branches do, so the result is the split set of the
corresponding unrooted tree — though summing is inexact in floating point, so lengths
recovered this way may differ in the last digit from the same tree written unrooted.
Lengths are `NaN` where the tree records none.

Splits are returned in a canonical order that depends only on the split set, so two trees
sharing a taxon index iterate their common splits identically.
"""
splits(tree; trivial::Bool = false) = splits(tree, taxonindex(tree); trivial)

function splits(tree, index::TaxonIndex{L}; trivial::Bool = false) where {L}
    ntaxa = length(index)
    masks = BitVector[]
    lengths = Dict{BitVector,Float64}()

    function record!(below, len)
        mask = _canonical(below)

        # A branch spanning every taxon separates nothing.
        count(mask) == 0 && return nothing
        !trivial && istrivial(mask, ntaxa) && return nothing

        if haskey(lengths, mask)
            lengths[mask] += len
        else
            lengths[mask] = len
            push!(masks, mask)
        end
        return nothing
    end

    _walksplits!(record!, tree, index, ntaxa)

    # Sorting makes the representation a function of the split set alone: two trees with
    # the same splits iterate them in the same order whatever their shapes.
    sort!(masks)

    return Splits{L}(index, masks, lengths)
end

taxonindex(s::Splits) = s.index

Base.length(s::Splits) = length(s.masks)
Base.isempty(s::Splits) = isempty(s.masks)
Base.eltype(::Type{<:Splits}) = BitVector
Base.iterate(s::Splits, state...) = iterate(s.masks, state...)
Base.pairs(s::Splits) = (mask => s.lengths[mask] for mask in s.masks)

Base.getindex(s::Splits, mask::AbstractVector{Bool}) = s.lengths[mask]
Base.get(s::Splits, mask::AbstractVector{Bool}, default) = get(s.lengths, mask, default)
Base.haskey(s::Splits, mask::AbstractVector{Bool}) = haskey(s.lengths, mask)
Base.in(mask::AbstractVector{Bool}, s::Splits) = haskey(s, mask)

function Base.show(io::IO, s::Splits)
    print(io, "Splits(", length(s), " over ", length(s.index), " taxa)")
end

# Masks from different taxon orderings index different taxa, so comparing them would
# silently answer a question nobody asked.
function _checksharedindex(a::Splits, b::Splits)
    a.index == b.index && return nothing
    throw(ArgumentError(
        "splits were built against different taxon orderings and cannot be compared; " *
        "build both with the index returned by taxonindex(t1, t2)"
    ))
end

for op in (:intersect, :union, :setdiff, :symdiff)
    @eval function Base.$op(a::Splits, b::Splits)
        _checksharedindex(a, b)
        return $op(a.masks, b.masks)
    end
end

"""
    incidencematrix(s::Splits) -> BitMatrix

A taxa × splits matrix whose `[i, j]` entry says whether taxon `i` lies on the marked side
of split `j`.

Rows follow the canonical taxon order, so the first row is always all `false`: splits are
oriented away from the first taxon.
"""
function incidencematrix(s::Splits)
    M = falses(length(s.index), length(s))
    for (j, mask) in enumerate(s.masks)
        M[:, j] = mask
    end
    return M
end
