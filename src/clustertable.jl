"""
    ClusterTable

The clusters of a tree, encoded so that membership can be tested in constant time.

Rooting a tree at one taxon and numbering its leaves depth-first makes every cluster — the
leaf set below some branch — a *contiguous* run of numbers. A cluster is then a pair of
endpoints rather than a set, and comparing two trees costs time proportional to their size
instead of to the product of their split counts.

Each cluster occupies one of two rows, `L` or `R`, and a lookup checks both. That two rows
always suffice follows from clusters being nested or disjoint: a cluster sharing its upper
endpoint with its parent cannot also share its lower one unless the parent has a single
child, so one row is always free.

Rooting at the taxon numbered first matches the orientation splits already use, so a
cluster is exactly a split read as the side excluding that taxon.

Day, W.H.E. (1985). *Optimal algorithms for comparing trees with labeled leaves.* Journal
of Classification 2(1): 7–28.
"""
struct ClusterTable
    code::Vector{Int32}     # code[p] is the leaf number of the taxon at index position p
    rowL::Vector{Int32}
    rowR::Vector{Int32}
    nclusters::Int
end

nclusters(table::ClusterTable) = table.nclusters

"""
    isclust(table, L, R) -> Bool

Whether the leaves numbered `L` through `R` form a cluster of the encoded tree.
"""
function isclust(table::ClusterTable, L::Int32, R::Int32)
    table.rowL[L] == L && table.rowR[L] == R && return true
    return table.rowL[R] == L && table.rowR[R] == R
end

"""
A tree flattened into parallel arrays: parents, children in compressed form, and the label
each leaf carries.

Working from indices rather than nodes lets the traversals run without recursion, so a tree
of any depth is safe, and keeps them free of pointer chasing. Labels are collected by the
same walk that records the structure, so reading a tree costs one traversal rather than one
per thing wanted from it.
"""
struct FlatTree{L}
    parent::Vector{Int32}       # 0 at the tree's own root
    childstart::Vector{Int32}   # childlist[childstart[i]:childstart[i+1]-1] are i's children
    childlist::Vector{Int32}
    leafat::Vector{Int32}       # node index -> position among the leaves; 0 if internal
    labels::Vector{L}           # by leaf position, in the order the walk reached them
end

nleaves(flat::FlatTree) = length(flat.labels)

"""
    flatten(tree) -> FlatTree

Read `tree` into flat arrays in a single walk.

The stack is typed by the node, so nodes of one tree must share a type — which is what any
concrete tree representation gives.
"""
function flatten(tree::T) where {T}
    parent = Int32[]
    leafat = Int32[]
    labels = Vector{_labeltype(T)}()

    nodes = T[tree]
    ups = Int32[0]

    # Depth-first and iterative: a node is discovered before its descendants, so a child
    # always lands at a higher position than its parent.
    while !isempty(nodes)
        node = pop!(nodes)
        up = pop!(ups)
        push!(parent, up)
        here = Int32(length(parent))

        kids = AbstractTrees.children(node)
        if isempty(kids)
            push!(labels, taxonlabel(node))
            push!(leafat, Int32(length(labels)))
        else
            push!(leafat, Int32(0))
            for kid in kids
                push!(nodes, kid)
                push!(ups, here)
            end
        end
    end

    nnodes = length(parent)

    # Compressed children, built by counting rather than by growing a vector per node.
    counts = zeros(Int32, nnodes)
    for up in parent
        up == 0 || (counts[up] += Int32(1))
    end
    childstart = Vector{Int32}(undef, nnodes + 1)
    childstart[1] = 1
    for i in 1:nnodes
        childstart[i + 1] = childstart[i] + counts[i]
    end
    next = copy(childstart)
    childlist = Vector{Int32}(undef, max(nnodes - 1, 0))
    for (i, up) in enumerate(parent)
        up == 0 && continue
        childlist[next[up]] = Int32(i)
        next[up] += Int32(1)
    end

    return FlatTree(parent, childstart, childlist, leafat, labels)
end

# The label type is whatever `taxonlabel` yields for this node type; asking the compiler
# keeps the label vector concrete without a leaf in hand.
_labeltype(::Type{T}) where {T} = Base.promote_op(taxonlabel, T)

# A node's neighbours are its children together with its parent. Rooting the tree elsewhere
# only changes which of them is "up", so every traversal here treats edges as undirected.
@inline function _neighbours(flat::FlatTree, i::Int32)
    kids = view(flat.childlist, flat.childstart[i]:(flat.childstart[i + 1] - 1))
    return kids, flat.parent[i]
end

"""
    taxonindex(f1::FlatTree, f2::FlatTree) -> TaxonIndex

The shared taxon ordering of two already-flattened trees, validated as the tree form is.

Flattening collects labels, so this needs no further walk.
"""
function taxonindex(f1::FlatTree, f2::FlatTree)
    _rejectrepeats(f1.labels)
    _rejectrepeats(f2.labels)
    return _sharedindex(f1.labels, f2.labels)
end

taxonindex(flat::FlatTree) = (_rejectrepeats(flat.labels); TaxonIndex(flat.labels))

"""
Map each leaf of a flattened tree to its position in the taxon ordering.
"""
_taxonpositions(flat::FlatTree, index::TaxonIndex) =
    Int32[index[label] for label in flat.labels]

"""
Number one tree's taxa by the order its own walk reached them, and number a second tree's
taxa to match.

Cluster comparison needs only that both trees agree on which taxon is which, never that the
numbering means anything outside the comparison. Sorting the labels — which is what a
[`TaxonIndex`](@ref) does, so that split masks are reproducible — would cost more than the
comparison itself, so it is skipped here. Both trees are still required to span the same
taxa, and a repeated label is still rejected.
"""
function _matchedpositions(f1::FlatTree{L}, f2::FlatTree{L}) where {L}
    n = length(f1.labels)
    length(f2.labels) == n || _reportdifferingtaxa(f1.labels, f2.labels)

    where = Dict{L,Int32}()
    sizehint!(where, n)
    for (i, label) in enumerate(f1.labels)
        haskey(where, label) && _rejectrepeats(f1.labels)
        where[label] = Int32(i)
    end

    second = Vector{Int32}(undef, n)
    taken = falses(n)
    for (i, label) in enumerate(f2.labels)
        at = get(where, label, Int32(0))
        at == 0 && _reportdifferingtaxa(f1.labels, f2.labels)
        taken[at] && _rejectrepeats(f2.labels)
        taken[at] = true
        second[i] = at
    end

    # The first tree's leaves are numbered by their own order, so their mapping is trivial.
    return Int32.(1:n), second
end

"""
Walk the tree depth-first from the leaf carrying taxon position `target` (default `1`),
returning the nodes in discovery order along with each node's parent under that rooting.
"""
function _rootedorder(flat::FlatTree, taxonpos::Vector{Int32}, target::Int32 = Int32(1))
    nnodes = Int32(length(flat.parent))

    start = Int32(0)
    for i in Int32(1):nnodes
        k = flat.leafat[i]
        if k != 0 && taxonpos[k] == target
            start = i
            break
        end
    end
    start == 0 && throw(ArgumentError(
        "tree does not contain the taxon the leaf numbering starts from"
    ))

    order = Vector{Int32}(undef, nnodes)
    up = Vector{Int32}(undef, nnodes)
    stack = Vector{Int32}(undef, nnodes)

    stack[1] = start
    up[start] = 0
    nstack = 1
    seen = 0

    while nstack > 0
        node = stack[nstack]
        nstack -= 1
        seen += 1
        order[seen] = node

        kids, parent = _neighbours(flat, node)
        from = up[node]
        for kid in kids
            kid == from && continue
            nstack += 1
            stack[nstack] = kid
            up[kid] = node
        end
        if parent != 0 && parent != from
            nstack += 1
            stack[nstack] = parent
            up[parent] = node
        end
    end

    return order, up
end

"""
Number the leaves in the order the walk reaches them, which is what makes each cluster a
contiguous run. The result is keyed by taxon position so that a second tree can be numbered
to match.
"""
function _daynumbers(flat::FlatTree, order::Vector{Int32}, taxonpos::Vector{Int32},
                     ntaxa::Int)
    code = Vector{Int32}(undef, ntaxa)
    seen = Int32(0)
    for node in order
        k = flat.leafat[node]
        k == 0 && continue
        seen += Int32(1)
        code[taxonpos[k]] = seen
    end
    return code
end

"""
Find the leaf interval, leaf count and number of branches below every node, working up from
the deepest.

Reversing the discovery order reaches a node only once all its descendants are known, so a
single backward pass suffices and no recursion is needed.
"""
function _intervals(flat::FlatTree, order::Vector{Int32}, up::Vector{Int32},
                    taxonpos::Vector{Int32}, code::Vector{Int32})
    nnodes = length(order)
    lo = Vector{Int32}(undef, nnodes)
    hi = Vector{Int32}(undef, nnodes)
    size = Vector{Int32}(undef, nnodes)
    below = Vector{Int32}(undef, nnodes)

    for pos in nnodes:-1:1
        node = order[pos]
        k = flat.leafat[node]

        if k != 0
            c = code[taxonpos[k]]
            lo[node] = c
            hi[node] = c
            size[node] = Int32(1)
            below[node] = Int32(0)
            continue
        end

        L = typemax(Int32)
        R = Int32(0)
        N = Int32(0)
        down = Int32(0)
        kids, parent = _neighbours(flat, node)
        from = up[node]
        for kid in kids
            kid == from && continue
            L = min(L, lo[kid])
            R = max(R, hi[kid])
            N += size[kid]
            down += Int32(1)
        end
        if parent != 0 && parent != from
            L = min(L, lo[parent])
            R = max(R, hi[parent])
            N += size[parent]
            down += Int32(1)
        end

        lo[node] = L
        hi[node] = R
        size[node] = N
        below[node] = down
    end

    return lo, hi, size, below
end

# A cluster worth recording is one that some tree on these taxa could lack: it must hold at
# least two taxa and leave at least two out, counting the taxon the numbering starts from.
# A node with a single branch below it describes the same cluster as that branch — which is
# what the root of a rooted tree becomes once the tree is rooted at a taxon instead — so it
# is passed over rather than counted twice.
@inline function _isinformative(N::Int32, below::Int32, ntaxa::Int)
    return below >= 2 && N >= 2 && N <= ntaxa - 2
end

"""
    clustertable(flat, index) -> ClusterTable
    clustertable(tree, index) -> ClusterTable

Encode a tree's clusters for constant-time membership tests.
"""
function clustertable(flat::FlatTree, taxonpos::Vector{Int32}, ntaxa::Int)
    order, up = _rootedorder(flat, taxonpos)
    code = _daynumbers(flat, order, taxonpos, ntaxa)
    lo, hi, size, below = _intervals(flat, order, up, taxonpos, code)

    rowL = zeros(Int32, ntaxa)
    rowR = zeros(Int32, ntaxa)
    total = 0

    for node in order
        flat.leafat[node] == 0 || continue
        parent = up[node]
        parent == 0 && continue

        _isinformative(size[node], below[node], ntaxa) || continue
        total += 1

        # Which endpoint the parent shares decides the row, leaving the other one free.
        L, R = lo[node], hi[node]
        row = hi[parent] == R ? L : R
        rowL[row] = L
        rowR[row] = R
    end

    return ClusterTable(code, rowL, rowR, total)
end

clustertable(flat::FlatTree, index::TaxonIndex) =
    clustertable(flat, _taxonpositions(flat, index), length(index))

clustertable(tree, index::TaxonIndex) = clustertable(flatten(tree), index)

"""
    sharedclusters(table, flat, index) -> (shared, nclusters)

How many of a tree's clusters the encoded tree also has, and how many it has in total.

A cluster can only match if its leaves are contiguous under the encoded tree's numbering,
which is a comparison of two integers; only then is the table consulted.
"""
function sharedclusters(table::ClusterTable, flat::FlatTree, taxonpos::Vector{Int32},
                        ntaxa::Int)
    order, up = _rootedorder(flat, taxonpos)
    lo, hi, size, below = _intervals(flat, order, up, taxonpos, table.code)

    shared = 0
    total = 0

    for node in order
        flat.leafat[node] == 0 || continue
        up[node] == 0 && continue

        N = size[node]
        _isinformative(N, below[node], ntaxa) || continue
        total += 1

        L, R = lo[node], hi[node]
        N == R - L + Int32(1) || continue
        isclust(table, L, R) && (shared += 1)
    end

    return shared, total
end

sharedclusters(table::ClusterTable, flat::FlatTree, index::TaxonIndex) =
    sharedclusters(table, flat, _taxonpositions(flat, index), length(index))

sharedclusters(table::ClusterTable, tree, index::TaxonIndex) =
    sharedclusters(table, flatten(tree), index)
