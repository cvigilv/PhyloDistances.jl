# The generalized Robinson-Foulds reduction: every metric in the split-matching family
# (Nye similarity, Jaccard-Robinson-Foulds, and — in later chunks — matching split
# distance and the information-theoretic metrics) differs from the others only in how a
# pair of splits is scored. This file supplies the shared matching step.

"""
    splitmatching(scorer, splits1::Splits, splits2::Splits; maximize::Bool = true)

The optimal matching between `splits1` and `splits2`, scored pairwise by
`scorer(mask1, mask2, ntaxa) -> Real` and reduced to a total by an assignment that
maximizes (or, with `maximize = false`, minimizes) the sum of matched scores.

`splits1` and `splits2` must share a taxon index (see [`splits`](@ref)); every split in
the smaller set is matched to a distinct split in the larger one, and any excess splits in
the larger set go unmatched. This is the shared substrate of the generalized Robinson-Foulds
family: [`NyeSimilarity`](@ref) and [`JaccardRobinsonFoulds`](@ref) are instances of it,
differing only in `scorer`, and a caller may supply its own `scorer` to define a new one.

Returns a `NamedTuple` `(score, matching)`: `score` is the total of the matched pairs'
scores, and `matching[i]` is the position in `splits2` matched to `splits1`'s `i`-th
split, or `0` if it goes unmatched.
"""
function splitmatching(scorer, splits1::Splits, splits2::Splits; maximize::Bool = true)
    _checksharedindex(splits1, splits2)
    ntaxa = length(taxonindex(splits1))
    n1, n2 = length(splits1), length(splits2)

    pairscore = Matrix{Float64}(undef, n1, n2)
    for j in 1:n2, i in 1:n1
        pairscore[i, j] = scorer(splits1.masks[i], splits2.masks[j], ntaxa)
    end

    return _matchtotal(pairscore; maximize)
end

# Shared tail of the matching pipeline: given an already-built pairwise score matrix,
# find the optimal matching and reduce it to a total. Split out so a metric that can score
# every pair faster than one scorer call at a time — see `_jaccardscorematrix` — can build
# `pairscore` itself and skip straight to this step.
function _matchtotal(pairscore::AbstractMatrix{Float64}; maximize::Bool)
    maximize && (pairscore .*= -1)
    rowmatch, _, total = _hungarian(pairscore)
    score = maximize ? -total : total
    return (score = score, matching = rowmatch)
end
