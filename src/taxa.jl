"""
    taxonlabel(leaf)

The label identifying `leaf` among the taxa of its tree.

Leaf labels are not part of the AbstractTrees.jl interface, so this is the hook a node
type must implement to be usable here. Labels must be sortable and hashable; they are
compared for equality across trees, so they carry the correspondence between two trees'
taxa.
"""
function taxonlabel(leaf)
    throw(ArgumentError(
        "PhyloDistances.taxonlabel is not defined for $(typeof(leaf)); " *
        "define it to return the taxon label of a leaf node"
    ))
end

taxonlabel(leaf::NewickTree.Node) = NewickTree.name(leaf)

"""
    TaxonIndex{L}

A canonical ordering of taxon labels, assigning each a position `1:length(index)`.

Splits, path matrices, and every other array a metric builds are indexed by these
positions, so two trees compared under one `TaxonIndex` agree on what row `i` means.
Labels are stored sorted, making the ordering a function of the taxon set alone rather
than of the order leaves happened to appear in either tree.

Construct with [`taxonindex`](@ref); index it with a label to recover a position.
"""
struct TaxonIndex{L}
    labels::Vector{L}
    positions::Dict{L,Int}

    TaxonIndex{L}(labels, positions) where {L} = new{L}(labels, positions)
end

function TaxonIndex(labels::AbstractVector)
    sorted = sort(collect(labels))
    return TaxonIndex{eltype(sorted)}(
        sorted, Dict(label => i for (i, label) in enumerate(sorted))
    )
end

"""
    taxa(index) -> Vector

The taxon labels of `index`, in canonical (sorted) order. Position `i` of the result is
the taxon that index `i` refers to.
"""
taxa(index::TaxonIndex) = index.labels

Base.length(index::TaxonIndex) = length(index.labels)
Base.getindex(index::TaxonIndex, label) = index.positions[label]
Base.haskey(index::TaxonIndex, label) = haskey(index.positions, label)
Base.:(==)(a::TaxonIndex, b::TaxonIndex) = a.labels == b.labels

function Base.show(io::IO, index::TaxonIndex)
    print(io, "TaxonIndex(", length(index), " taxa: ")
    join(io, repr.(first(index.labels, 4)), ", ")
    length(index) > 4 && print(io, ", …")
    print(io, ")")
end

"""
    taxonlabels(tree) -> Vector

The labels of `tree`'s leaves, in traversal order.

Throws if a label occurs more than once. A repeated label would make the correspondence
between two trees' taxa ambiguous, and every downstream computation would silently use
whichever leaf happened to be visited last.
"""
function taxonlabels(tree)
    labels = [taxonlabel(leaf) for leaf in AbstractTrees.Leaves(tree)]

    counts = Dict{eltype(labels),Int}()
    for label in labels
        counts[label] = get(counts, label, 0) + 1
    end
    repeated = sort!([label for (label, n) in counts if n > 1])

    isempty(repeated) || throw(ArgumentError(
        "tree has repeated taxon labels: $(join(repr.(repeated), ", ")); " *
        "taxon labels must be unique within a tree"
    ))
    return labels
end

"""
    taxonindex(tree) -> TaxonIndex
    taxonindex(t1, t2) -> TaxonIndex

The canonical taxon ordering of a tree, or the shared ordering of two trees.

The two-tree form additionally requires that both trees span the same taxon set, and
throws naming the labels that differ. Metrics use it to obtain their indexing, so calling
it is what validates that a comparison is meaningful.
"""
taxonindex(tree) = TaxonIndex(taxonlabels(tree))

function taxonindex(t1, t2)
    l1, l2 = taxonlabels(t1), taxonlabels(t2)
    s1, s2 = Set(l1), Set(l2)
    if s1 != s2
        only1 = sort!(collect(setdiff(s1, s2)))
        only2 = sort!(collect(setdiff(s2, s1)))
        throw(ArgumentError(
            "trees span different taxa: " *
            "only in the first tree: $(isempty(only1) ? "none" : join(repr.(only1), ", "))" *
            "; only in the second tree: " *
            "$(isempty(only2) ? "none" : join(repr.(only2), ", "))"
        ))
    end
    return TaxonIndex(l1)
end

"""
    isrooted(tree) -> Bool

Whether `tree` is rooted, read from the number of children at its root: two means rooted,
three or more means the root is an arbitrary starting point for writing the tree down and
carries no evolutionary claim. This is the usual reading of Newick, which has no explicit
marker for rooting.

Throws for a root with fewer than two children, where the question has no answer.
"""
function isrooted(tree)
    n = length(AbstractTrees.children(tree))
    n == 2 && return true
    n >= 3 && return false
    throw(ArgumentError(
        "cannot tell whether a tree is rooted from a root with $n " *
        (n == 1 ? "child" : "children") * "; a tree needs at least three taxa"
    ))
end

# Reconciles the rooting of the inputs against what the comparison is defined on, so a
# metric's own implementation may assume it receives trees of the rooting it declares.
function _checkrooting(comparison, t1, t2)
    wanted = requiresrooted(comparison)
    _checkonerooting(comparison, t1, wanted, "first")
    _checkonerooting(comparison, t2, wanted, "second")
    return nothing
end

function _checkonerooting(comparison, tree, wanted::Bool, position::AbstractString)
    rooted = isrooted(tree)
    rooted == wanted && return nothing

    name = nameof(typeof(comparison))
    if rooted
        # The root's two child branches induce the same bipartition, so discarding the
        # root position leaves a well-defined unrooted tree.
        @warn "$name is defined on unrooted trees but the $position tree is rooted; " *
              "its root position is ignored and its splits are read as unrooted " *
              "bipartitions."
    else
        throw(ArgumentError(
            "$name is defined on rooted trees but the $position tree is unrooted. " *
            "There is no canonical way to root it — midpoint and outgroup rooting give " *
            "different answers — so root it explicitly before comparing."
        ))
    end
    return nothing
end
