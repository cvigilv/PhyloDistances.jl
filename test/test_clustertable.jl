using PhyloDistances
using PhyloDistances: FlatTree, ClusterTable, clustertable, flatten, isclust, nclusters,
    sharedclusters
using Random: Xoshiro
using Test

# The split-set formulation of Robinson-Foulds, kept as the definition the interval
# encoding has to reproduce.
splitrf(t1, t2) = (i = taxonindex(t1, t2); length(symdiff(splits(t1, i), splits(t2, i))))

"""Collapse an unrooted tree's degree-3 root to degree 2, making it rooted."""
function asrooted(tree)
    kids = collect(PhyloDistances.NewickTree.children(tree))
    length(kids) < 3 && return tree
    root = PhyloDistances.NewickTree.Node(UInt16(9990), PhyloDistances.NewickTree.NewickData())
    rest = PhyloDistances.NewickTree.Node(
        UInt16(9991), PhyloDistances.NewickTree.NewickData(1.0, 0.0, "")
    )
    push!(root, kids[1])
    for kid in kids[2:end]
        push!(rest, kid)
    end
    push!(root, rest)
    return root
end

@testset "flatten" begin
    flat = flatten(readnw("(A,B,(C,D));"))

    @test PhyloDistances.nleaves(flat) == 4
    @test sort(flat.labels) == ["A", "B", "C", "D"]

    # One entry per node: four leaves, the root, and the internal node above C and D.
    @test length(flat.parent) == 6
    @test count(iszero, flat.parent) == 1        # only the tree's own root has no parent
    @test count(!iszero, flat.leafat) == 4

    @testset "children are reachable through the compressed lists" begin
        for i in 1:length(flat.parent)
            kids = flat.childlist[flat.childstart[i]:(flat.childstart[i + 1] - 1)]
            @test all(flat.parent[k] == i for k in kids)
        end
        @test sum(flat.childstart[end] - flat.childstart[1]) == length(flat.parent) - 1
    end
end

@testset "cluster counts match the split count" begin
    rng = Xoshiro(3)
    for n in 4:14
        tree = randomtree(rng, n)
        flat = flatten(tree)
        table = clustertable(flat, taxonindex(flat))

        @test nclusters(table) == length(splits(tree)) == n - 3
    end
end

@testset "isclust finds exactly the tree's clusters" begin
    # ((C,D) is a cluster; leaves are numbered from A, which roots the walk.
    tree = readnw("(A,B,(C,D));")
    flat = flatten(tree)
    index = taxonindex(flat)
    table = clustertable(flat, index)

    @test nclusters(table) == 1

    # Exactly one interval of length two is a cluster, and it is the only one recorded.
    found = [(L, R) for L in Int32(1):Int32(4), R in Int32(1):Int32(4) if
             L <= R && isclust(table, L, R)]
    @test length(found) == 1
end

@testset "a rooted tree does not double-count its root" begin
    # Rooting at a taxon leaves the original root with one branch below it, describing the
    # same cluster as that branch. Counting both would inflate every distance.
    rooted = readnw("((A,B),(C,D));")
    unrooted = readnw("(A,B,(C,D));")

    for tree in (rooted, unrooted)
        flat = flatten(tree)
        @test nclusters(clustertable(flat, taxonindex(flat))) == 1
    end

    @test RobinsonFoulds()(rooted, rooted) == 0
    @test splitrf(rooted, unrooted) == 0
end

@testset "agrees with the split-set definition" begin
    @testset "unrooted trees" begin
        rng = Xoshiro(11)
        for _ in 1:200
            n = rand(rng, 4:30)
            a = randomtree(rng, n)
            b = perturb(rng, a, rand(rng, 0:(2n)))
            @test RobinsonFoulds()(a, b) == splitrf(a, b)
        end
    end

    @testset "rooted trees, where the root becomes a pass-through node" begin
        rng = Xoshiro(12)
        for _ in 1:200
            n = rand(rng, 4:30)
            base = randomtree(rng, n)
            # Perturb before rooting: rooting re-parents the tree's own children.
            moved = perturb(rng, base, rand(rng, 0:(2n)))
            a, b = asrooted(base), asrooted(moved)
            @test RobinsonFoulds()(a, b) == splitrf(a, b)
        end
    end

    @testset "trees with polytomies carry fewer clusters" begin
        star = readnw("(" * join(("T$i" for i in 1:8), ",") * ");")
        flat = flatten(star)
        @test nclusters(clustertable(flat, taxonindex(flat))) == 0
        @test RobinsonFoulds()(star, star) == 0

        partly = readnw("(T1,T2,T3,T4,(T5,T6,T7,T8));")
        pflat = flatten(partly)
        @test nclusters(clustertable(pflat, taxonindex(pflat))) == length(splits(partly))
    end
end

@testset "sharedclusters" begin
    a = readnw("(T1,T2,((T3,T4),T5));")
    b = readnw("(T1,T2,((T3,T5),T4));")
    fa, fb = flatten(a), flatten(b)
    index = taxonindex(fa, fb)
    table = clustertable(fa, index)

    shared, total = sharedclusters(table, fb, index)
    @test total == length(splits(b))
    @test shared == length(intersect(splits(a, index), splits(b, index)))

    @testset "a tree against itself shares everything" begin
        same, n = sharedclusters(table, fa, index)
        @test same == n == nclusters(table)
    end
end

@testset "deep trees do not overflow the stack" begin
    # The traversals are iterative, so depth costs memory rather than stack frames.
    rng = Xoshiro(21)
    tall = randomtree(rng, 3000)
    @test RobinsonFoulds()(tall, tall) == 0
    @test RobinsonFoulds()(tall, perturb(rng, tall, 1)) == 2
end
