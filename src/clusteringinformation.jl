# Mutual clustering information and its complementary distance. Both compare split sets
# through the same optimal assignment, with mutual information as the pair score.

# Sum the binary entropy carried by a tree's non-trivial splits. The taxon count comes from
# the shared index rather than from each mask separately.
function _treeclusteringentropy(s::Splits, logcounts::Vector{Float64})
    n = length(taxonindex(s))
    return sum(mask -> _clusteringentropy(logcounts, count(mask), n), s; init = 0.0)
end

function _treeclusteringentropy(s::Splits)
    return _treeclusteringentropy(s, _log2counts(length(taxonindex(s))))
end

# log2 of every integer contingency count. Mutual information only takes logarithms of
# counts between zero and the taxon count, so one table replaces millions of transcendental
# calls when scoring all pairs of splits. The zero entry is finite because a zero-count cell
# contributes zero and never needs log2(0).
function _log2counts(n::Integer)
    logs = Vector{Float64}(undef, n + 1)
    logs[1] = 0.0
    for k in 1:n
        logs[k + 1] = log2(k)
    end
    return logs
end

# Binary entropy expressed through integer logarithms. Keeping tree entropy and exact-match
# scores on this same path preserves exact self-comparisons.
function _clusteringentropy(logcounts::Vector{Float64}, k::Integer, n::Integer)
    return (
        n * logcounts[n + 1] - k * logcounts[k + 1] -
            (n - k) * logcounts[n - k + 1]
    ) / n
end

# Build the mutual-information score matrix directly. Marginal split sizes are invariant
# across a row or column, so their counts and logarithmic terms are computed once. The
# remaining logarithms depend only on integer contingency counts and come from `logcounts`.
function _mutualinformationscorematrix(
    masks1::Vector{BitVector}, masks2::Vector{BitVector}, ntaxa::Integer,
    logcounts::Vector{Float64}
)
    w1, w2 = _packmasks(masks1), _packmasks(masks2)
    na = [count(mask) for mask in masks1]
    nb = [count(mask) for mask in masks2]
    pairscore = Matrix{Float64}(undef, length(masks1), length(masks2))
    nwords = size(w1, 1)

    logn = logcounts[ntaxa + 1]
    offseta = [logn - logcounts[k + 1] for k in na]
    offsetA = [logn - logcounts[ntaxa - k + 1] for k in na]
    logb = [logcounts[k + 1] for k in nb]
    logB = [logcounts[ntaxa - k + 1] for k in nb]

    # The matrix axes also index the matching marginal and offset vectors, and every
    # contingency count lies in 0:ntaxa, so each indexed collection is covered here.
    @inbounds for j in axes(pairscore, 2)
        nbj = nb[j]
        for i in axes(pairscore, 1)
            nai = na[i]
            nA = ntaxa - nai
            aandb = _countand(w1, w2, i, j, nwords)

            # Equal and complementary masks induce the same partition. Use the entropy
            # path so self-information remains exact, including after exact-match removal
            # when this matrix helper is called directly elsewhere.
            if aandb == nai == nbj || (aandb == 0 && nai + nbj == ntaxa)
                pairscore[i, j] = _clusteringentropy(logcounts, nai, ntaxa)
                continue
            end

            aandB = nai - aandb
            Aandb = nbj - aandb
            AandB = nA - Aandb

            # Four equal cells describe independent balanced partitions, whose mutual
            # information is exactly zero. This also avoids cancellation noise.
            if aandb == aandB == Aandb == AandB
                pairscore[i, j] = 0.0
                continue
            end

            total =
                aandb * (logcounts[aandb + 1] + offseta[i] - logb[j]) +
                aandB * (logcounts[aandB + 1] + offseta[i] - logB[j]) +
                Aandb * (logcounts[Aandb + 1] + offsetA[i] - logb[j]) +
                AandB * (logcounts[AandB + 1] + offsetA[i] - logB[j])
            pairscore[i, j] = total / ntaxa
        end
    end
    return pairscore
end

function _mutualinformationscorematrix(
    masks1::Vector{BitVector}, masks2::Vector{BitVector}, ntaxa::Integer
)
    return _mutualinformationscorematrix(masks1, masks2, ntaxa, _log2counts(ntaxa))
end

_mutualinformationscorematrix(s1::Splits, s2::Splits, ntaxa::Integer) =
    _mutualinformationscorematrix(s1.masks, s2.masks, ntaxa)

# Some optimum always matches an exact split pair. If either split's previous partner was
# unmatched, MI(A, X) <= H(A) makes the exchange non-decreasing. If both had other partners,
# the triangle inequality for variation of information shows that matching the exact pair
# and the two displaced partners cannot lower the total MI. Remove exact pairs before
# constructing the dense score matrix and assignment problem, then add their entropy back.
function _unmatchedsplits(
    s1::Splits, s2::Splits, ntaxa::Integer, logcounts::Vector{Float64}
)
    unmatched1 = Int[]
    unmatched2 = Int[]
    exactscore = 0.0

    for (i, mask) in pairs(s1.masks)
        if haskey(s2, mask)
            exactscore += _clusteringentropy(logcounts, count(mask), ntaxa)
        else
            push!(unmatched1, i)
        end
    end
    for (i, mask) in pairs(s2.masks)
        haskey(s1, mask) || push!(unmatched2, i)
    end

    return unmatched1, unmatched2, exactscore
end

_unmatchedsplits(s1::Splits, s2::Splits, ntaxa::Integer) =
    _unmatchedsplits(s1, s2, ntaxa, _log2counts(ntaxa))

function _mutualclusteringinfo(
    s1::Splits, s2::Splits, ntaxa::Integer, logcounts::Vector{Float64}
)
    unmatched1, unmatched2, exactscore = _unmatchedsplits(s1, s2, ntaxa, logcounts)

    # If either remainder is empty, every split on the smaller side was matched exactly;
    # unmatched excess splits contribute no similarity.
    (isempty(unmatched1) || isempty(unmatched2)) && return exactscore

    masks1 = length(unmatched1) == length(s1) ? s1.masks : s1.masks[unmatched1]
    masks2 = length(unmatched2) == length(s2) ? s2.masks : s2.masks[unmatched2]
    pairscore = _mutualinformationscorematrix(masks1, masks2, ntaxa, logcounts)
    score = exactscore + _matchtotal(pairscore; maximize = true).score
    return iszero(score) ? 0.0 : score
end

_mutualclusteringinfo(s1::Splits, s2::Splits, ntaxa::Integer) =
    _mutualclusteringinfo(s1, s2, ntaxa, _log2counts(ntaxa))

"""
    MutualClusteringInfo(; convention = :treedist, normalize = false)

The mutual clustering information similarity of Smith (2020), in bits. Each split is
viewed as a two-class partition of the taxa. The similarity is the largest total
[`mutualinformation`](@ref) over a one-to-one matching between the two trees' non-trivial
splits. Excess splits in either tree remain unmatched and contribute zero.

The unnormalized self-similarity of a tree is the sum of
[`clusteringentropy`](@ref) over its splits. With `normalize = true`, the result is divided
by the mean of the two trees' summed split entropies, following TreeDist. Identical trees
therefore have normalized similarity 1. A function such as `max` or `min` instead combines
the two per-tree entropy totals, and a number is used directly as the divisor.

Both conventions currently compute the same definition. The comparison removes identical
split pairs, then constructs the remaining pair scores and solves a dense linear assignment
problem.

Smith, M.R. (2020). *Information theoretic Generalized Robinson-Foulds metrics for
comparing phylogenetic trees.* Bioinformatics 36(20): 5007–5013.

See [`ClusteringInfoDistance`](@ref) for the complementary distance and
[`splitmatching`](@ref) for the matching framework.
"""
struct MutualClusteringInfo{C <: Convention, N} <: TreeSimilarity
    convention::C
    normalize::N
end

MutualClusteringInfo(; convention = TreeDistConvention(), normalize = false) =
    MutualClusteringInfo(Convention(convention), normalize)

function _compare(::MutualClusteringInfo, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    s1, s2 = splits(t1, index), splits(t2, index)
    return _mutualclusteringinfo(s1, s2, length(index))
end

normalizerinfo(::MutualClusteringInfo, ::Convention, tree) =
    _treeclusteringentropy(splits(tree, taxonindex(tree)))

function normalizer(comparison::MutualClusteringInfo, convention::Convention, t1, t2)
    return (
        normalizerinfo(comparison, convention, t1) +
            normalizerinfo(comparison, convention, t2)
    ) / 2
end

"""
    ClusteringInfoDistance(; convention = :treedist, normalize = false)

The clustering information distance of Smith (2020), in bits. It converts
[`MutualClusteringInfo`](@ref) to a distance by subtracting twice the optimal mutual
clustering information from the summed clustering entropy of both trees:

```math
CID(T_1, T_2) = H(T_1) + H(T_2) - 2 MCI(T_1, T_2).
```

Here `H(T)` is the sum of [`clusteringentropy`](@ref) over the tree's non-trivial splits.
This is the variation of information between optimally matched split partitions, with an
unmatched split contributing its full entropy. The result is symmetric, vanishes for
identical trees, and satisfies the triangle inequality.

With `normalize = true`, the distance is divided by `H(T_1) + H(T_2)`, following TreeDist.
A function such as `max` or `min` instead combines the two entropy totals, and a number is
used directly as the divisor. Both conventions currently compute the same definition.

Smith, M.R. (2020). *Information theoretic Generalized Robinson-Foulds metrics for
comparing phylogenetic trees.* Bioinformatics 36(20): 5007–5013.

Meilă, M. (2007). *Comparing clusterings: an information based distance.* Journal of
Multivariate Analysis 98(5): 873–895.
"""
struct ClusteringInfoDistance{C <: Convention, N} <: TreeMetric
    convention::C
    normalize::N
end

ClusteringInfoDistance(; convention = TreeDistConvention(), normalize = false) =
    ClusteringInfoDistance(Convention(convention), normalize)

function _compare(::ClusteringInfoDistance, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    s1, s2 = splits(t1, index), splits(t2, index)
    ntaxa = length(index)
    logcounts = _log2counts(ntaxa)
    entropy =
        _treeclusteringentropy(s1, logcounts) + _treeclusteringentropy(s2, logcounts)
    return entropy - 2 * _mutualclusteringinfo(s1, s2, ntaxa, logcounts)
end

normalizerinfo(::ClusteringInfoDistance, ::Convention, tree) =
    _treeclusteringentropy(splits(tree, taxonindex(tree)))
