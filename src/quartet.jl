"""
    QuartetDistance(; convention = :treedist, normalize = false)

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

# Complexity

`O(n²)` to tabulate the leaf-to-leaf path lengths of each tree, then `O(n⁴)` to enumerate
the `binomial(n, 4)` quartets, each resolved in constant time. The quartet count grows
steeply: 210 subsets at 10 taxa, 64 million at 200, 41 billion at 1000. This is the exact
distance by direct enumeration, and it is the right choice only at modest taxon counts.

Estabrook, G.F., McMorris, F.R. and Meacham, C.A. (1985). *Comparison of undirected
phylogenetic trees based on subtrees of four evolutionary units.* Systematic Zoology
34(2): 193–200.
"""
struct QuartetDistance{C<:Convention,N} <: TreeMetric
    convention::C
    normalize::N
end

QuartetDistance(; convention = TreeDistConvention(), normalize = false) =
    QuartetDistance(Convention(convention), normalize)

function _compare(::QuartetDistance, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
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
