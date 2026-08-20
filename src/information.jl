"""
    log2rooted(n) -> Float64

log2 of the number of rooted binary trees on `n` labeled leaves, `(2n - 3)!!`.

Computed as a running sum of `log2(2j - 3)` rather than the double factorial itself, so it
stays finite where `(2n - 3)!!` would overflow any integer or `Float64` representation —
already true past `n ≈ 300`. By convention `log2rooted(0) == log2rooted(1) == 0`: there is
exactly one (trivial) rooted tree on 0 or 1 leaves.
"""
function log2rooted(n::Integer)
    n < 0 && throw(ArgumentError("log2rooted is undefined for negative n, got $n"))
    n <= 1 && return 0.0
    total = 0.0
    for j in 2:n
        total += log2(2j - 3)
    end
    return total
end

"""
    log2unrooted(n) -> Float64

log2 of the number of unrooted binary trees on `n` labeled leaves, `(2n - 5)!!` for
`n ≥ 3`.

Related to [`log2rooted`](@ref) by `log2unrooted(n) = log2rooted(n) - log2(2n - 3)`, since
unrooting a rooted tree on `n` leaves collapses `2n - 3` equally likely root positions onto
one unrooted topology. By convention `log2unrooted(n) == 0` for `n ≤ 2`.
"""
function log2unrooted(n::Integer)
    n < 0 && throw(ArgumentError("log2unrooted is undefined for negative n, got $n"))
    n <= 2 && return 0.0
    return log2rooted(n) - log2(2n - 3)
end

"""
    splitinfo(k, n) -> Float64
    splitinfo(mask::AbstractVector{Bool}) -> Float64

The phylogenetic information content of a split separating `k` of `n` taxa from the rest:
`-log2` of the probability that a uniformly chosen unrooted binary tree on `n` taxa
contains it, `log2unrooted(n) - log2rooted(k) - log2rooted(n - k)`.

A split shared by every tree on the same taxa — one side holding a single taxon, i.e. a
trivial split — carries exactly zero information, for any `n`; passing a `mask` gives `k`
and `n` as `count(mask)` and `length(mask)`.

Thorley, J.L., Wilkinson, M. and Charleston, M. (1998). *The information content of
consensus trees.* In Advances in Data Science and Classification, 91–98. Springer.
"""
function splitinfo(k::Integer, n::Integer)
    (0 <= k <= n) || throw(ArgumentError("splitinfo needs 0 <= k <= n, got k=$k, n=$n"))
    return log2unrooted(n) - log2rooted(k) - log2rooted(n - k)
end

splitinfo(mask::AbstractVector{Bool}) = splitinfo(count(mask), length(mask))

"""
    clusteringentropy(k, n) -> Float64
    clusteringentropy(mask::AbstractVector{Bool}) -> Float64

The Shannon entropy, in bits, of a split treated as a two-class partition of `n` taxa into
groups of size `k` and `n - k`: `-p*log2(p) - (1-p)*log2(1-p)` with `p = k/n`. This is the
number of bits needed to encode which side of the split a uniformly drawn taxon falls on;
it is `0` where `p ∈ {0, 1}` (every taxon on one side) and maximal, `1` bit, at `p = 1/2`.

Passing a `mask` gives `k` and `n` as `count(mask)` and `length(mask)`.

Meilă, M. (2007). *Comparing clusterings—an information based distance.* Journal of
Multivariate Analysis 98(5): 873–895.
"""
function clusteringentropy(k::Integer, n::Integer)
    (0 <= k <= n) || throw(ArgumentError(
        "clusteringentropy needs 0 <= k <= n, got k=$k, n=$n"
    ))
    (n == 0 || k == 0 || k == n) && return 0.0
    p = k / n
    return -(p * log2(p) + (1 - p) * log2(1 - p))
end

clusteringentropy(mask::AbstractVector{Bool}) = clusteringentropy(count(mask), length(mask))

"""
    mutualinformation(mask1::AbstractVector{Bool}, mask2::AbstractVector{Bool}) -> Float64

The mutual information, in bits, between the two-class partitions two splits induce on a
shared taxon set: how many bits knowing which side of `mask1` a randomly drawn taxon falls
on reveals about which side of `mask2` it falls on.

Computed from the 2×2 contingency table of the four overlap counts (both marked, `mask1`
marked only, `mask2` marked only, neither marked) as `sum(x * log2(x * n / (row * col)))`
over the nonzero cells, then divided by `n` to convert the extensive sum into a per-taxon
quantity on the same scale as [`clusteringentropy`](@ref). A split's mutual information
with itself equals its own entropy.

Vinh, N.X., Epps, J. and Bailey, J. (2010). *Information theoretic measures for
clusterings comparison: variants, properties, normalization and correction for chance.*
Journal of Machine Learning Research 11: 2837–2854.
"""
function mutualinformation(mask1::AbstractVector{Bool}, mask2::AbstractVector{Bool})
    n = length(mask1)
    length(mask2) == n || throw(DimensionMismatch(
        "mutualinformation needs both masks over the same taxa: got lengths " *
        "$(length(mask1)) and $(length(mask2))"
    ))
    n == 0 && return 0.0

    na, nb = count(mask1), count(mask2)
    nA, nB = n - na, n - nb

    aandb = 0
    for i in eachindex(mask1, mask2)
        aandb += mask1[i] & mask2[i]
    end
    aandB = na - aandb
    Aandb = nb - aandb
    AandB = nA - Aandb

    total = 0.0
    for (x, row, col) in ((aandb, na, nb), (aandB, na, nB), (Aandb, nA, nb), (AandB, nA, nB))
        x == 0 && continue
        total += x * log2(x * n / (row * col))
    end
    return total / n
end

"""
    jointentropy(mask1::AbstractVector{Bool}, mask2::AbstractVector{Bool}) -> Float64

The joint entropy, in bits, of the two-class partitions two splits induce on a shared taxon
set: `clusteringentropy(mask1) + clusteringentropy(mask2) - mutualinformation(mask1, mask2)`.
"""
function jointentropy(mask1::AbstractVector{Bool}, mask2::AbstractVector{Bool})
    return clusteringentropy(mask1) + clusteringentropy(mask2) -
           mutualinformation(mask1, mask2)
end
