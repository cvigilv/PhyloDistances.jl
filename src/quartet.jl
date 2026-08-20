"""
    QuartetDistance(; convention = :treedist, normalize = false, algorithm = :fast)

The quartet distance: the number of four-taxon subsets the two trees resolve differently.

Any four taxa are related by one of three unrooted topologies, `ab|cd`, `ac|bd` or `ad|bc`,
unless the tree leaves them unresolved. Two trees are compared by tallying the four-taxon
subsets on which they disagree, which counts both a subset the trees resolve in *conflicting*
ways and one that only a single tree resolves at all. A fully resolved tree is therefore at
the maximal distance `binomial(n, 4)` from the star tree on the same taxa. Some authors
instead count only direct conflicts, treating an unresolved quartet as agreeing with
everything; the two definitions coincide on binary trees, where nothing is unresolved.

Branch lengths are ignored, and the trees are read as unrooted.

With `normalize = true` the result is divided by `binomial(n, 4)`, the number of four-taxon
subsets, so the value is the fraction of quartets on which the trees disagree.

Both conventions compute the same value. TreeDist has no quartet distance — it belongs to
the companion R package [Quartet](https://github.com/ms609/Quartet) — so there is no
reference formulation to diverge from.

# Algorithm

`algorithm = :fast` (the default) counts concordant quartets without enumerating them, in
`O(n³)`; see [`_fastconcordantcount`](@ref) for the method. It requires that at least one
of the two trees be fully resolved (binary); if both carry a polytomy, it falls back to
`:naive` and warns, since exactness cannot be guaranteed without the general (and
unimplemented) two-polytomy case.

`algorithm = :naive` always uses direct enumeration: `O(n²)` to tabulate the leaf-to-leaf
path lengths of each tree, then `O(n⁴)` to enumerate the `binomial(n, 4)` quartets, each
resolved in constant time. It remains the correctness oracle the fast path is tested
against. The quartet count grows steeply: 210 subsets at 10 taxa, 64 million at 200,
41 billion at 1000 — `:naive` is the right choice only at modest taxon counts.

Estabrook, G.F., McMorris, F.R. and Meacham, C.A. (1985). *Comparison of undirected
phylogenetic trees based on subtrees of four evolutionary units.* Systematic Zoology
34(2): 193–200.
"""
struct QuartetDistance{C<:Convention,N} <: TreeMetric
    convention::C
    normalize::N
    algorithm::Symbol

    function QuartetDistance{C,N}(convention, normalize, algorithm) where {C,N}
        algorithm in (:fast, :naive) || throw(ArgumentError(
            "unknown algorithm $(repr(algorithm)); expected :fast or :naive"
        ))
        return new{C,N}(convention, normalize, algorithm)
    end
end

QuartetDistance(convention::C, normalize::N, algorithm::Symbol = :fast) where {C,N} =
    QuartetDistance{C,N}(convention, normalize, algorithm)

QuartetDistance(; convention = TreeDistConvention(), normalize = false, algorithm = :fast) =
    QuartetDistance(Convention(convention), normalize, algorithm)

function _compare(m::QuartetDistance, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    length(index) < 4 && return 0
    m.algorithm === :naive && return _naivequartetdistance(t1, t2, index)

    fast = _fastquartetdistance(t1, t2, index)
    fast === nothing || return fast

    @warn "QuartetDistance's :fast algorithm needs at least one input tree to be fully " *
          "resolved; both trees have a polytomy, so falling back to the O(n⁴) :naive " *
          "enumeration, which will be slow on large trees."
    return _naivequartetdistance(t1, t2, index)
end

function _naivequartetdistance(t1, t2, index::TaxonIndex)
    n = length(index)
    d1 = _topologicalpaths(t1, index)
    d2 = _topologicalpaths(t2, index)

    differing = 0
    for a in 1:(n - 3), b in (a + 1):(n - 2), c in (b + 1):(n - 1), d in (c + 1):n
        q1 = _quartettopology(
            d1[a, b] + d1[c, d], d1[a, c] + d1[b, d], d1[a, d] + d1[b, c]
        )
        q2 = _quartettopology(
            d2[a, b] + d2[c, d], d2[a, c] + d2[b, d], d2[a, d] + d2[b, c]
        )
        differing += q1 != q2
    end
    return differing
end

"""
    _fastquartetdistance(t1, t2, index::TaxonIndex) -> Union{Int,Nothing}

The quartet distance computed without enumerating quartets, or `nothing` if neither tree
is fully resolved (see [`_fastconcordantcount`](@ref) for why one binary tree is required).

For binary trees every quartet is resolved, so the distance is `binomial(n, 4)` minus the
number resolved identically by both trees, `_fastconcordantcount`. The same holds whenever
at least one tree is binary: a quartet unresolved by a polytomous tree can only be
concordant if the binary tree is *also* unresolved there, which never happens.
"""
function _fastquartetdistance(t1, t2, index::TaxonIndex)
    n = length(index)
    f1, f2 = flatten(t1), flatten(t2)
    pos1 = Int32[index[label] for label in f1.labels]
    pos2 = Int32[index[label] for label in f2.labels]

    _isbinaryflat(f1, pos1) || _isbinaryflat(f2, pos2) || return nothing

    concordant = _fastconcordantcount(f1, f2, pos1, pos2, n)
    return binomial(n, 4) - concordant
end

# Whether every internal node has exactly two branches below it, under an arbitrary
# rooting (branch count below a node does not depend on which taxon the tree is rooted
# at, so checking one rooting settles it for all of them).
function _isbinaryflat(flat::FlatTree, positions::Vector{Int32})
    n = length(flat.labels)
    n < 3 && return true
    order, up = _rootedorder(flat, positions, Int32(1))
    code = _daynumbers(flat, order, positions, n)
    _, _, _, below = _intervals(flat, order, up, positions, code)
    return all(b -> b <= 2, below)
end

"""
    _fastconcordantcount(f1, f2, pos1, pos2, n) -> Int

The number of four-taxon subsets `t1` and `t2` resolve identically, computed in `O(n³)`.

# Method

A quartet `{a,b,c,d}` resolved as `ab|cd` is, from `a`'s point of view, a *rooted triple*:
rooting the tree at `a` makes `b` and `c,d`'s most recent common ancestor a strict
descendant of `a,b,c,d`'s overall ancestor, i.e. `b` and one of `{c,d}` are not the closest
pair — concretely, `{b,c,d}` resolves as `(cd)b`. Every unrooted quartet is a rooted triple
under each of its four members in turn, so summing rooted-triple agreement between `t1` and
`t2` over all `n` choices of root and dividing by four gives the quartet agreement count.

Fix a root `w` and the remaining `n - 1` leaves. A pair `x, y` has a most recent common
ancestor `v` in the tree rooted at `w`; write `clade(v)` for the leaves below `v`. The
triple `{x, y, z}` resolves as `(xy)z` for every `z` outside `clade(v)`, and this is the
*only* way `{x, y}` can be the close pair — so summing, over every unordered pair `x, y`,
the count of `z` for which `x, y` are the close pair in *both* trees gives the number of
concordant triples rooted at `w`, and each such `z` lies outside `clade₁(v₁) ∪ clade₂(v₂)`
for `v₁ = mrca₁(x,y)`, `v₂ = mrca₂(x,y)`.

Rooting `t1` at `w` and numbering the other `n - 1` leaves by that walk's discovery order
makes every clade a contiguous run of numbers (as in [`clustertable`](@ref)); the same is
done for `t2`, under an independent numbering. For each branch point `v₁` of `t1`, a dense
array over `t2`'s numbering records which leaves lie in `clade₁(v₁)`, and its prefix sum
turns "how many of `clade₁(v₁)`'s leaves lie in `clade₂(v₂)`" into one O(1) lookup — used
once per pair `x, y` whose most recent common ancestor in `t1` is `v₁`, of which there are
`O(n)` in total. Building the array costs `O(n)` and is paid once per branch point, so one
root costs `O(n²)`; summing over the `n` roots gives the `O(n³)` bound. A branch point
serving too few pairs to repay that `O(n)` counts each of its pairs directly instead, in
`|clade₁(v₁)|` steps, which is bounded by the same `O(n²)` per root.
"""
function _fastconcordantcount(
    f1::FlatTree, f2::FlatTree, pos1::Vector{Int32}, pos2::Vector{Int32}, n::Int
)
    # One row of padding. `lca2` is read and written at scattered `[x, y]`, so a leading
    # dimension that is a power of two sends whole families of those accesses to the same
    # cache set: at 1024 taxa the unpadded matrix runs the whole count two to four times
    # slower than the padded one. The extra row is never indexed.
    lca2 = Matrix{UInt64}(undef, n + 1, n)
    sigma = Vector{Int32}(undef, n)         # sigma[t1 position] == t2 position, same taxon
    tau = Vector{Int32}(undef, n)           # tau[t2 position] == t1 position, same taxon
    prefix = Vector{Int32}(undef, n + 1)
    ranges = Tuple{Int32,Int32}[]

    total = 0
    for w in Int32(1):Int32(n)
        order1, up1 = _rootedorder(f1, pos1, w)
        code1 = _daynumbers(f1, order1, pos1, n)
        lo1, hi1, _, _ = _intervals(f1, order1, up1, pos1, code1)

        order2, up2 = _rootedorder(f2, pos2, w)
        code2 = _daynumbers(f2, order2, pos2, n)
        lo2, hi2, _, _ = _intervals(f2, order2, up2, pos2, code2)
        _fillcrosspairs!(lca2, ranges, f2, order2, up2, lo2, hi2)

        for taxon in 1:n
            sigma[code1[taxon]] = code2[taxon]
            tau[code2[taxon]] = code1[taxon]
        end

        total += _crosspaircontribution!(
            prefix, ranges, f1, order1, up1, lo1, hi1, sigma, tau, lca2, n
        )
    end

    quotient, remainder = divrem(total, 4)
    remainder == 0 || throw(ErrorException(
        "PhyloDistances internal error: concordant-quartet accumulator " *
        "$total is not divisible by 4"
    ))
    return quotient
end

# A clade's interval packed into a single word, its two endpoints in the halves of a
# `UInt64`. The two are always written and read together, at indices that jump around a
# matrix of `n²` entries, so carrying them in one word halves the memory traffic that
# dominates both loops below.
_packinterval(lo::Int32, hi::Int32) = (UInt64(lo) << 32) | UInt64(hi)
_unpackinterval(word::UInt64) = (Int32(word >> 32), Int32(word & 0xffffffff))

# Fills lca[x, y] with the interval spanned by the most recent common ancestor of the
# leaves numbered x and y, for every pair under this rooting. A node's downward neighbours
# (its children, plus its `FlatTree` parent if that direction has become "down" under this
# rooting) partition the leaves below it into contiguous, disjoint runs; every pair drawn
# from two different runs has its most recent common ancestor here, and every pair has its
# most recent common ancestor at exactly one node, so this reaches each of the O(n²) pairs
# once.
function _fillcrosspairs!(
    lca::Matrix{UInt64}, ranges::Vector{Tuple{Int32,Int32}}, flat::FlatTree,
    order::Vector{Int32}, up::Vector{Int32}, lo::Vector{Int32}, hi::Vector{Int32}
)
    for node in order
        kids, parent = _neighbours(flat, node)
        from = up[node]
        empty!(ranges)
        for kid in kids
            kid == from || push!(ranges, (lo[kid], hi[kid]))
        end
        parent != 0 && parent != from && push!(ranges, (lo[parent], hi[parent]))

        k = length(ranges)
        k < 2 && continue
        word = _packinterval(lo[node], hi[node])
        for i in 1:(k - 1), j in (i + 1):k
            loi, hii = ranges[i]
            loj, hij = ranges[j]
            for x in loi:hii, y in loj:hij
                lca[x, y] = word
                lca[y, x] = word
            end
        end
    end
    return nothing
end

# How many of the leaves numbered `L1:H1` in one tree's ordering fall within `l2lo:l2hi`
# in the other's, read through `sigma`, which translates between the two numberings.
function _countshared(
    sigma::Vector{Int32}, L1::Int32, H1::Int32, l2lo::Int32, l2hi::Int32
)
    shared = Int32(0)
    for q in L1:H1
        s = sigma[q]
        shared += (l2lo <= s <= l2hi) ? Int32(1) : Int32(0)
    end
    return shared
end

# For the tree rooted at `w`, adds up (n - 1 - |clade1(mrca1(x,y)) ∪ clade2(mrca2(x,y))|)
# over every pair x, y of the other n - 1 leaves, grouped by their most recent common
# ancestor in `f1` (found the same way `_fillcrosspairs!` finds it in `f2`) so that the
# `t2`-indexed membership array for `clade1` can be built once per branch point rather than
# once per pair.
function _crosspaircontribution!(
    prefix::Vector{Int32}, ranges::Vector{Tuple{Int32,Int32}}, f1::FlatTree,
    order1::Vector{Int32}, up1::Vector{Int32}, lo1::Vector{Int32}, hi1::Vector{Int32},
    sigma::Vector{Int32}, tau::Vector{Int32}, lca2::Matrix{UInt64}, n::Int
)
    contribution = 0
    for node in order1
        kids, parent = _neighbours(f1, node)
        from = up1[node]
        empty!(ranges)
        for kid in kids
            kid == from || push!(ranges, (lo1[kid], hi1[kid]))
        end
        parent != 0 && parent != from && push!(ranges, (lo1[parent], hi1[parent]))

        k = length(ranges)
        k < 2 && continue

        L1, H1 = lo1[node], hi1[node]
        clade1size = H1 - L1 + 1

        npairs = 0
        for i in 1:(k - 1), j in (i + 1):k
            npairs += (ranges[i][2] - ranges[i][1] + 1) * (ranges[j][2] - ranges[j][1] + 1)
        end

        # Tabulating costs n whatever the node serves, while answering one pair directly
        # costs this clade's own size: a cherry deep in the tree would scan every taxon to
        # serve its single pair. Tabulate only where that is the cheaper of the two, which
        # keeps the O(n³) bound and drops most of the constant.
        tabulate = npairs * clade1size > n
        if tabulate
            # prefix[p + 1] counts, among the first p leaves of *t2*'s numbering, how many
            # also lie in clade1(node) — found via `tau`, which reads their t1 position.
            acc = Int32(0)
            prefix[1] = Int32(0)
            for p in 1:n
                acc += (L1 <= tau[p] <= H1) ? Int32(1) : Int32(0)
                prefix[p + 1] = acc
            end
        end

        for i in 1:(k - 1), j in (i + 1):k
            loi, hii = ranges[i]
            loj, hij = ranges[j]
            for x in loi:hii
                x2 = sigma[x]
                for y in loj:hij
                    l2lo, l2hi = _unpackinterval(lca2[x2, sigma[y]])
                    clade2size = l2hi - l2lo + 1
                    intersection = tabulate ? prefix[l2hi + 1] - prefix[l2lo] :
                                              _countshared(sigma, L1, H1, l2lo, l2hi)
                    unionsize = clade1size + clade2size - intersection
                    contribution += n - 1 - unionsize
                end
            end
        end
    end
    return contribution
end

# Both trees span the same taxa, so each defines the same set of quartets; `max` and `min`
# combine to that same count rather than to a different scale.
normalizerinfo(::QuartetDistance, ::Convention, tree) =
    binomial(length(taxonlabels(tree)), 4)

normalizer(::QuartetDistance, ::Convention, t1, t2) =
    binomial(length(taxonindex(t1, t2)), 4)

Distances.result_type(m::QuartetDistance, ::Type, ::Type) =
    isnormalized(m) ? Float64 : Int

"""
    _quartettopology(abcd, acbd, adbc) -> Int

Which of the three pairings of four taxa a tree resolves, given the three sums of
opposite-pair path lengths: `1` for `ab|cd`, `2` for `ac|bd`, `3` for `ad|bc`, and `0`
where the tree resolves none of them.

By the four-point condition the sum belonging to the induced topology is the smallest, and
the other two are equal, exceeding it by twice the length of the path separating the two
pairs. All three coincide exactly when no such path exists, which is the unresolved case.
"""
function _quartettopology(abcd::Integer, acbd::Integer, adbc::Integer)
    abcd < acbd && abcd < adbc && return 1
    acbd < abcd && acbd < adbc && return 2
    adbc < abcd && adbc < acbd && return 3
    return 0
end

"""
    _topologicalpaths(tree, index::TaxonIndex) -> Matrix{Int}

The number of edges between each pair of leaves, indexed by taxon position.

Every pair of leaves meets at exactly one node — their most recent common ancestor — so
recording a pair's distance where its two subtrees join reaches each pair once and builds
the whole matrix in `O(n²)`.

A node with a single branch below it, such as the root of a rooted tree, subdivides an edge
and so inflates the distances across it. Quartet topologies are unaffected: subdividing
leaves the path lengths a tree metric with positive edge weights, which is all the
four-point condition requires.
"""
function _topologicalpaths(tree, index::TaxonIndex)
    n = length(index)
    distances = zeros(Int, n, n)
    depths = zeros(Int, n)
    _pathsbelow!(distances, depths, tree, index, 0)
    return distances
end

# Returns the taxon positions of the leaves below `node`, filling in the distance between
# every pair of leaves whose most recent common ancestor is `node`.
function _pathsbelow!(distances, depths, node, index::TaxonIndex, depth::Int)
    kids = AbstractTrees.children(node)

    if isempty(kids)
        position = index[taxonlabel(node)]
        depths[position] = depth
        return [position]
    end

    below = Int[]
    for child in kids
        subtree = _pathsbelow!(distances, depths, child, index, depth + 1)
        for i in subtree, j in below
            distances[i, j] = depths[i] + depths[j] - 2 * depth
            distances[j, i] = distances[i, j]
        end
        append!(below, subtree)
    end
    return below
end
