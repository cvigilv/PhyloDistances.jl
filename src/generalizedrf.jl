# Nye similarity and the Jaccard-Robinson-Foulds distance: the split-matching family's two
# best-known instances, both scored by the same Jaccard-index pair scorer and differing
# only in the exponent applied to it and in whether conflicting splits may be matched.

# Population count of `a .& b` without materializing the intermediate BitVector that
# broadcasting `.&` would allocate. A split-matching metric calls this O(n²) times — once
# per pair of splits — so the allocation avoided here is the single largest cost in
# building its score matrix. Safe with no bounds handling: BitVector's unused trailing
# bits are always zero, and `a`/`b` share their taxon count by construction (both are
# built from the same TaxonIndex), so their chunk vectors are the same length.
function _countand(a::BitVector, b::BitVector)
    s = 0
    for i in eachindex(a.chunks, b.chunks)
        s += count_ones(a.chunks[i] & b.chunks[i])
    end
    return s
end

# The pair scorer TreeDist calls jaccard_similarity: for splits a|A of tree 1 and b|B of
# tree 2, four Jaccard indices are available depending on which side of each split is
# paired with which. min(J(a,b), J(A,B)) reads the pairing "a corresponds to b" (and so A
# to B); min(J(a,B), J(A,b)) reads the opposite pairing. Taking the larger of the two picks
# whichever orientation places the splits nearer each other. Only an identical pair reaches
# a Jaccard index of 1 in both halves of that pairing at once — a merely nested pair (one
# properly contained in the other) always leaves at least one of the four ratios below 1 —
# so raising the result to `k` sharpens it toward an indicator of exact agreement as `k`
# grows: a mismatched pair's score falls monotonically toward zero, with no special case
# needed for `k = Inf`.
#
# Split into a count-based core and a mask-based wrapper so that scoring every pair of two
# trees' splits (`_jaccardscorematrix`) can hoist each split's marked-taxon count out of
# the O(n²) loop — it depends on the row or column alone, never on the pair — while a
# standalone call still takes masks directly.
function _jaccardscorecore(
    a_and_b::Integer, na::Integer, nb::Integer, ntaxa::Integer; k::Real, allowconflict::Bool
)
    nA, nB = ntaxa - na, ntaxa - nb
    a_and_B = na - a_and_b
    A_and_b = nb - a_and_b
    A_and_B = nB - a_and_B

    if !allowconflict
        # Compatible iff one split's marked (or unmarked) side nests entirely inside a
        # side of the other — the standard split-compatibility condition.
        compatible = a_and_b == na || a_and_B == na || A_and_b == nA || A_and_B == nA
        compatible || return 0.0
    end

    jaccard_ab = a_and_b / (ntaxa - A_and_B)
    jaccard_AB = A_and_B / (ntaxa - a_and_b)
    jaccard_aB = a_and_B / (ntaxa - A_and_b)
    jaccard_Ab = A_and_b / (ntaxa - a_and_B)

    similarity = max(min(jaccard_ab, jaccard_AB), min(jaccard_aB, jaccard_Ab))
    return similarity^k
end

function _jaccardscore(a::BitVector, b::BitVector, ntaxa::Integer; k::Real, allowconflict::Bool)
    return _jaccardscorecore(_countand(a, b), count(a), count(b), ntaxa; k, allowconflict)
end

# The score matrix a split-matching metric optimizes over, built directly rather than
# through `splitmatching`'s one-scorer-call-per-cell interface: every split's
# marked-taxon count is computed once here instead of once per pair scored against it.
function _jaccardscorematrix(
    s1::Splits, s2::Splits, ntaxa::Integer; k::Real, allowconflict::Bool
)
    n1, n2 = length(s1), length(s2)
    na = [count(m) for m in s1.masks]
    nb = [count(m) for m in s2.masks]

    pairscore = Matrix{Float64}(undef, n1, n2)
    for j in 1:n2
        b = s2.masks[j]
        nbj = nb[j]
        for i in 1:n1
            pairscore[i, j] =
                _jaccardscorecore(_countand(s1.masks[i], b), na[i], nbj, ntaxa; k, allowconflict)
        end
    end
    return pairscore
end

"""
    NyeSimilarity(; convention = :treedist, normalize = false)

The similarity of Nye _et al._ (2006): the optimal matching between two trees' splits,
scored by the size of the largest split consistent with both, normalized against the
Jaccard index.

Equivalent to [`JaccardRobinsonFoulds`](@ref) with `k = 1` and `allowconflict = true`,
computed as a similarity (largest on identical trees) rather than a distance.

With `normalize = true` the result is divided by the mean of the two trees' split counts,
which reaches exactly 1 when the trees are identical.

Nye, T.M.W., Liò, P. and Gascuel, O. (2006). *A novel algorithm and web-based tool for
comparing two alternative phylogenetic trees.* Bioinformatics 22(1): 117–119.

See [`splitmatching`](@ref) for the matching framework this reduces to, and
[`JaccardRobinsonFoulds`](@ref) for the distance form with a tunable exponent.
"""
struct NyeSimilarity{C<:Convention,N} <: TreeSimilarity
    convention::C
    normalize::N
end

NyeSimilarity(; convention = TreeDistConvention(), normalize = false) =
    NyeSimilarity(Convention(convention), normalize)

function _compare(::NyeSimilarity, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    s1, s2 = splits(t1, index), splits(t2, index)
    ntaxa = length(index)
    pairscore = _jaccardscorematrix(s1, s2, ntaxa; k = 1, allowconflict = true)
    return _matchtotal(pairscore; maximize = true).score
end

normalizerinfo(::NyeSimilarity, ::Convention, tree) = length(splits(tree, taxonindex(tree)))

function normalizer(comparison::NyeSimilarity, convention::Convention, t1, t2)
    return (normalizerinfo(comparison, convention, t1) +
            normalizerinfo(comparison, convention, t2)) / 2
end

"""
    JaccardRobinsonFoulds(; k = 1, allowconflict = true, convention = :treedist, normalize = false)

The Jaccard-Robinson-Foulds distance: the optimal matching between two trees' splits,
scored by the largest Jaccard index consistent with both and raised to the exponent `k`,
converted to a distance by doubling the matched total and subtracting it from the number
of splits the two trees carry between them.

Splits that conflict — cannot both appear on one tree — may still be matched by default;
set `allowconflict = false` to score every conflicting pair zero instead, so the matching
prefers to leave them unpaired. As `k` grows a matched pair's score falls toward zero
unless the splits agree exactly, so the distance rises monotonically toward
[`RobinsonFoulds`](@ref) and never exceeds it. `k = 1` with `allowconflict = true` gives
the same scoring as [`NyeSimilarity`](@ref), computed here as a distance rather than a
similarity.

With `normalize = true` the result is divided by the total number of splits the two trees
carry between them, as for [`RobinsonFoulds`](@ref).

Böcker, S., Canzar, S. and Klau, G.W. (2013). *The generalized Robinson-Foulds metric.* In
Algorithms in Bioinformatics (WABI 2013), Lecture Notes in Computer Science 8126: 156–169.

See [`splitmatching`](@ref) for the underlying framework and [`NyeSimilarity`](@ref) for
the `k = 1`, `allowconflict = true` case as a similarity.
"""
struct JaccardRobinsonFoulds{C<:Convention,N,K<:Real} <: TreeMetric
    convention::C
    normalize::N
    k::K
    allowconflict::Bool
end

function JaccardRobinsonFoulds(;
    convention = TreeDistConvention(), normalize = false, k = 1, allowconflict = true
)
    return JaccardRobinsonFoulds(Convention(convention), normalize, k, allowconflict)
end

function _compare(m::JaccardRobinsonFoulds, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    s1, s2 = splits(t1, index), splits(t2, index)
    ntaxa = length(index)
    pairscore = _jaccardscorematrix(s1, s2, ntaxa; k = m.k, allowconflict = m.allowconflict)
    similarity = _matchtotal(pairscore; maximize = true).score
    return length(s1) + length(s2) - 2 * similarity
end

normalizerinfo(::JaccardRobinsonFoulds, ::Convention, tree) =
    length(splits(tree, taxonindex(tree)))
