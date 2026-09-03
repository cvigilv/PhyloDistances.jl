# Mutual clustering information and its complementary distance. Both compare split sets
# through the same optimal assignment, with mutual information as the pair score.

# Sum the binary entropy carried by a tree's non-trivial splits. The taxon count comes from
# the shared index rather than from each mask separately.
function _treeclusteringentropy(s::Splits)
    n = length(taxonindex(s))
    return sum(mask -> clusteringentropy(count(mask), n), s; init = 0.0)
end

# Build the mutual-information score matrix directly. Marginal split sizes are invariant
# across a row or column, so counting each once avoids repeating that work for every pair.
function _mutualinformationscorematrix(s1::Splits, s2::Splits, ntaxa::Integer)
    w1, w2 = _packmasks(s1.masks), _packmasks(s2.masks)
    na = [count(mask) for mask in s1.masks]
    nb = [count(mask) for mask in s2.masks]
    pairscore = Matrix{Float64}(undef, length(s1), length(s2))
    nwords = size(w1, 1)

    for j in axes(pairscore, 2)
        nbj = nb[j]
        for i in axes(pairscore, 1)
            pairscore[i, j] = _mutualinformation(
                _countand(w1, w2, i, j, nwords), na[i], nbj, ntaxa
            )
        end
    end
    return pairscore
end

function _mutualclusteringinfo(s1::Splits, s2::Splits, ntaxa::Integer)
    pairscore = _mutualinformationscorematrix(s1, s2, ntaxa)
    score = _matchtotal(pairscore; maximize = true).score
    return iszero(score) ? 0.0 : score
end

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

Both conventions currently compute the same definition. The comparison constructs all
split-pair scores and solves a dense linear assignment problem.

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
    entropy = _treeclusteringentropy(s1) + _treeclusteringentropy(s2)
    return entropy - 2 * _mutualclusteringinfo(s1, s2, length(index))
end

normalizerinfo(::ClusteringInfoDistance, ::Convention, tree) =
    _treeclusteringentropy(splits(tree, taxonindex(tree)))
