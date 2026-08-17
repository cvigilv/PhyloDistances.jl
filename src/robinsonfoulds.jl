"""
    RobinsonFoulds(; convention = :treedist, normalize = false)

The Robinson-Foulds distance: the number of splits present in one tree but not the other.

Two trees are compared purely by which bipartitions of the taxa they contain, so branch
lengths are ignored and only the topology counts. The distance rises by two for each split
one tree has that the other lacks, and is zero exactly when the trees carry the same
splits. Trivial splits are excluded, since every tree on a taxon set has all of them.

With `normalize = true` the result is divided by the total number of splits the two trees
carry between them. That total is what the trees could differ by at most, so a normalized
distance of one means they share no split at all. For binary trees this equals `2(n - 3)`,
but the two part company as soon as either tree has a polytomy and so carries fewer splits.

Both conventions compute the same value; the definition is not in dispute.

Robinson, D.F. and Foulds, L.R. (1981). *Comparison of phylogenetic trees.* Mathematical
Biosciences 53(1–2): 131–147.

See [`WeightedRobinsonFoulds`](@ref) to take branch lengths into account.
"""
struct RobinsonFoulds{C<:Convention,N} <: TreeMetric
    convention::C
    normalize::N
end

RobinsonFoulds(; convention = TreeDistConvention(), normalize = false) =
    RobinsonFoulds(Convention(convention), normalize)

function _compare(::RobinsonFoulds, ::Convention, t1, t2)
    # One walk per tree yields structure and labels together, and the two are numbered
    # against each other rather than against a sorted ordering neither result depends on.
    f1, f2 = flatten(t1), flatten(t2)
    pos1, pos2 = _matchedpositions(f1, f2)
    n = length(f1.labels)

    table = clustertable(f1, pos1, n)
    shared, n2 = sharedclusters(table, f2, pos2, n)

    # Each tree's unmatched clusters are what separates them.
    return nclusters(table) + n2 - 2 * shared
end

function normalizerinfo(::RobinsonFoulds, ::Convention, tree)
    flat = flatten(tree)
    _rejectrepeats(flat.labels)
    n = length(flat.labels)
    return nclusters(clustertable(flat, Int32.(1:n), n))
end

Distances.result_type(m::RobinsonFoulds, ::Type, ::Type) =
    isnormalized(m) ? Float64 : Int

"""
    WeightedRobinsonFoulds(; convention = :treedist, normalize = false)

The weighted Robinson-Foulds distance: the total disagreement in branch lengths, summed
over every split either tree contains.

A split present in both trees contributes the difference between its two lengths, and one
present in only a single tree contributes its whole length, as though the other tree gave
it length zero. Trivial splits are included, so lengths on the branches leading to
individual taxa count like any other.

Unlike [`RobinsonFoulds`](@ref), two trees of identical topology are at a nonzero distance
when their branch lengths differ. Every branch must carry a length; a tree written without
them is rejected rather than silently compared as `NaN`.

With `normalize = true` the result is divided by the two trees' combined branch length,
which is the most they could disagree by.

Both conventions compute the same value. TreeDist has no branch-length-weighted
Robinson-Foulds, so there is no reference implementation to diverge from.

Robinson, D.F. and Foulds, L.R. (1979). *Comparison of weighted labelled trees.* In
Combinatorial Mathematics VI, Lecture Notes in Mathematics 748: 119–126.
"""
struct WeightedRobinsonFoulds{C<:Convention,N} <: TreeMetric
    convention::C
    normalize::N
end

WeightedRobinsonFoulds(; convention = TreeDistConvention(), normalize = false) =
    WeightedRobinsonFoulds(Convention(convention), normalize)

function _compare(::WeightedRobinsonFoulds, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    s1 = _weightedsplits(t1, index)
    s2 = _weightedsplits(t2, index)

    # A split missing from one tree is read as length zero there, so the sum runs over
    # every split either tree contains.
    return sum(union(s1, s2); init = 0.0) do mask
        abs(get(s1, mask, 0.0) - get(s2, mask, 0.0))
    end
end

normalizerinfo(::WeightedRobinsonFoulds, ::Convention, tree) =
    sum(last, pairs(_weightedsplits(tree, taxonindex(tree))); init = 0.0)

# Comparing lengths is meaningless when a tree records none, and NaN would carry that
# silently through every downstream sum.
function _weightedsplits(tree, index::TaxonIndex)
    s = splits(tree, index; trivial = true)
    any(isnan, last(pair) for pair in pairs(s)) && throw(ArgumentError(
        "a branch-length metric needs every branch to have a length, " *
        "and this tree leaves at least one unset"
    ))
    return s
end
