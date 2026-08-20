using PhyloDistances
using PhyloDistances: SplitKey, branchlength, incidencematrix, istrivial
using Random
using Test

"""The canonical mask over `index` marking exactly `labels`, oriented away from taxon 1."""
function mask(index, labels...)
    m = falses(length(index))
    for label in labels
        m[index[label]] = true
    end
    return first(m) ? .!m : m
end

"""A caterpillar tree on `n` taxa, unrooted: root of degree 3, all branches length 1."""
function caterpillar(n)
    nw = "T$n:1.0"
    for i in (n - 1):-1:3
        nw = "($nw,T$i:1.0):1.0"
    end
    return readnw("(T1:1.0,T2:1.0,$nw);")
end

@testset "branchlength" begin
    tree = readnw("((A:1,B:2):3,C:4);")
    @test branchlength(first(PhyloDistances.AbstractTrees.children(tree))) == 3.0

    # Newick records no length above the root, and none at all in a topology-only tree.
    @test isnan(branchlength(tree))
    @test isnan(branchlength(first(PhyloDistances.AbstractTrees.children(readnw("((A,B),C);")))))

    @test_throws "branchlength is not defined for Int64" branchlength(1)
end

@testset "istrivial" begin
    # Canonical masks never contain taxon 1, so one taxon marked or all-but-one marked
    # both describe a pendant branch.
    @test istrivial(BitVector([0, 1, 0, 0]), 4)
    @test istrivial(BitVector([0, 1, 1, 1]), 4)
    @test !istrivial(BitVector([0, 1, 1, 0]), 4)
end

@testset "canonical orientation" begin
    index = taxonindex(readnw("((A:1,B:2):3,(C:4,D:5):6);"))
    s = splits(readnw("((A:1,B:2):3,(C:4,D:5):6);"), index)

    # Splits are oriented away from the first taxon, so A is never a member.
    for m in s
        @test !m[index["A"]]
    end

    # {A,B} and {C,D} are the same bipartition and must be one split, not two.
    @test length(s) == 1
    @test only(s) == mask(index, "C", "D")
    @test only(s) == mask(index, "A", "B")
end

@testset "hand-computed split sets" begin
    @testset "four taxa, unrooted" begin
        tree = readnw("(T1:1.0,T2:1.0,(T3:1.0,T4:1.0):5.0);")
        index = taxonindex(tree)
        s = splits(tree, index)

        @test length(s) == 1
        @test only(s) == mask(index, "T3", "T4")
        @test s[mask(index, "T3", "T4")] == 5.0
    end

    @testset "five taxa, unrooted" begin
        tree = readnw("(T1:1.0,T2:1.0,((T3:1.0,T4:1.0):2.0,T5:1.0):3.0);")
        index = taxonindex(tree)
        s = splits(tree, index)

        @test length(s) == 2
        @test mask(index, "T3", "T4") in s
        @test mask(index, "T3", "T4", "T5") in s
        @test s[mask(index, "T3", "T4")] == 2.0
        @test s[mask(index, "T3", "T4", "T5")] == 3.0
    end

    @testset "three taxa have no informative splits" begin
        @test isempty(splits(readnw("(A:1,B:1,C:1);")))
    end
end

@testset "a rooted tree yields its unrooted split set" begin
    rooted = readnw("((T1:1.0,T2:1.0):3.0,(T3:1.0,T4:1.0):6.0);")
    unrooted = readnw("(T1:1.0,T2:1.0,(T3:1.0,T4:1.0):9.0);")

    index = taxonindex(rooted, unrooted)
    sr, su = splits(rooted, index), splits(unrooted, index)

    @test collect(sr) == collect(su)

    # The two root branches describe one branch of the unrooted tree, so their lengths add.
    @test sr[mask(index, "T3", "T4")] == 9.0
    @test sr[mask(index, "T3", "T4")] == su[mask(index, "T3", "T4")]

    @testset "summing lengths is inexact" begin
        # Recovering an unrooted branch by addition is a floating-point sum, so metrics
        # comparing lengths across rootings must not demand exact equality.
        r = readnw("((T1:1.0,T2:1.0):0.2,(T3:1.0,T4:1.0):0.4);")
        u = readnw("(T1:1.0,T2:1.0,(T3:1.0,T4:1.0):0.6);")
        i = taxonindex(r, u)
        m = mask(i, "T3", "T4")

        @test splits(r, i)[m] != splits(u, i)[m]
        @test splits(r, i)[m] ≈ splits(u, i)[m]
    end
end

@testset "iteration order depends only on the split set" begin
    # The same splits reached by different traversals must come out in the same order.
    a = readnw("(T1:1.0,T2:1.0,((T3:1.0,T4:1.0):1.0,T5:1.0):1.0);")
    b = readnw("((T3:1.0,T4:1.0):1.0,T5:1.0,(T1:1.0,T2:1.0):1.0);")
    index = taxonindex(a, b)

    @test collect(splits(a, index)) == collect(splits(b, index))
end

@testset "trivial splits" begin
    tree = readnw("(T1:1.0,T2:2.0,(T3:3.0,T4:4.0):5.0);")
    index = taxonindex(tree)

    @test length(splits(tree, index)) == 1
    @test length(splits(tree, index; trivial = true)) == 1 + 4

    withtrivial = splits(tree, index; trivial = true)
    @test withtrivial[mask(index, "T2")] == 2.0
    @test withtrivial[mask(index, "T3")] == 3.0

    # Taxon 1's pendant branch is the all-but-first mask, since orientation excludes it.
    @test withtrivial[mask(index, "T2", "T3", "T4")] == 1.0
end

@testset "split count invariants" begin
    @testset "a binary unrooted tree on n taxa has n-3 informative splits" begin
        for n in 4:9
            @test length(splits(caterpillar(n))) == n - 3
        end
    end

    @testset "including trivial splits adds one per taxon" begin
        for n in 4:9
            @test length(splits(caterpillar(n); trivial = true)) == (n - 3) + n
        end
    end

    @testset "a star tree has no informative splits" begin
        star = readnw("(" * join(("T$i:1.0" for i in 1:6), ",") * ");")
        @test isempty(splits(star))
        @test length(splits(star; trivial = true)) == 6
    end
end

@testset "splits depend on the taxon set, not on leaf order" begin
    a = readnw("(T1:1.0,T2:1.0,(T3:1.0,T4:1.0):5.0);")
    b = readnw("(T4:1.0,T3:1.0,(T2:1.0,T1:1.0):5.0);")
    index = taxonindex(a, b)

    # Same bipartition {T1,T2} = {T3,T4} written from either side.
    @test collect(splits(a, index)) == collect(splits(b, index))
end

@testset "orientation is idempotent" begin
    tree = readnw("(T1:1.0,T2:1.0,((T3:1.0,T4:1.0):2.0,T5:1.0):3.0);")
    index = taxonindex(tree)
    s = splits(tree, index)

    # Re-canonicalizing a stored mask must be a no-op.
    for m in s
        @test PhyloDistances._canonical(m) == m
    end
end

@testset "set operations" begin
    a = readnw("(T1:1.0,T2:1.0,((T3:1.0,T4:1.0):1.0,T5:1.0):1.0);")
    b = readnw("(T1:1.0,T2:1.0,((T4:1.0,T5:1.0):1.0,T3:1.0):1.0);")
    index = taxonindex(a, b)
    sa, sb = splits(a, index), splits(b, index)

    @test length(intersect(sa, sb)) == 1
    @test issetequal(intersect(sa, sb), [mask(index, "T3", "T4", "T5")])

    @test length(union(sa, sb)) == 3
    @test setdiff(sa, sb) == [mask(index, "T3", "T4")]
    @test setdiff(sb, sa) == [mask(index, "T4", "T5")]
    @test length(symdiff(sa, sb)) == 2

    @testset "identical trees" begin
        @test isempty(symdiff(sa, splits(a, index)))
        @test length(intersect(sa, splits(a, index))) == length(sa)
    end

    @testset "mismatched taxon orderings are rejected" begin
        other = readnw("(X1:1.0,X2:1.0,(X3:1.0,X4:1.0):1.0);")
        so = splits(other)

        @test_throws "different taxon orderings" symdiff(sa, so)
        @test_throws "taxonindex(t1, t2)" intersect(sa, so)
    end
end

@testset "lookup and iteration" begin
    tree = readnw("(T1:1.0,T2:1.0,(T3:1.0,T4:1.0):5.0);")
    index = taxonindex(tree)
    s = splits(tree, index)

    @test taxonindex(s) == index
    @test eltype(s) === BitVector
    @test collect(pairs(s)) == [mask(index, "T3", "T4") => 5.0]
    @test get(s, mask(index, "T2", "T3"), -1.0) == -1.0
    @test !haskey(s, mask(index, "T2", "T3"))
    @test_throws KeyError s[mask(index, "T2", "T3")]

    @test occursin("Splits(1 over 4 taxa)", sprint(show, s))
end

@testset "split identity is the underlying words" begin
    tree = readnw("(T1:1.0,T2:1.0,(T3:1.0,T4:1.0):5.0);")
    index = taxonindex(tree)
    s = splits(tree, index)
    m = mask(index, "T3", "T4")

    # A split is identified by its contents, not by the array type carrying them, so a
    # caller's `Vector{Bool}` finds the same split a `BitVector` does.
    @test s[collect(Bool, m)] == 5.0
    @test haskey(s, collect(Bool, m))

    # The unused bits of a mask's final word are zero, so a mask over more taxa can carry
    # exactly the same words while describing a different bipartition.
    longer = vcat(m, falses(1))
    @test m.chunks == longer.chunks
    @test !isequal(SplitKey(longer), SplitKey(m))
    @test !haskey(s, longer)
    @test_throws KeyError s[longer]
end

@testset "set operations match elementwise membership" begin
    # Splits are compared a word at a time; these sizes put the taxa one under, one over,
    # and exactly on a 64-bit word boundary, where zero padding would show up as a
    # disagreement with plain elementwise equality.
    rng = Xoshiro(20260820)
    for n in (8, 65, 128)
        a, b = randomtree(rng, n), randomtree(rng, n)
        index = taxonindex(a, b)
        sa, sb = splits(a, index), splits(b, index)
        ma, mb = collect(sa), collect(sb)

        inb(m) = any(==(m), mb)
        ina(m) = any(==(m), ma)

        @test intersect(sa, sb) == filter(inb, ma)
        @test setdiff(sa, sb) == filter(!inb, ma)
        @test union(sa, sb) == vcat(ma, filter(!ina, mb))
        @test symdiff(sa, sb) == vcat(filter(!inb, ma), filter(!ina, mb))
    end
end

@testset "incidencematrix" begin
    tree = readnw("(T1:1.0,T2:1.0,((T3:1.0,T4:1.0):2.0,T5:1.0):3.0);")
    index = taxonindex(tree)
    s = splits(tree, index)
    M = incidencematrix(s)

    @test size(M) == (5, 2)

    # Orientation excludes the first taxon from every split.
    @test !any(M[1, :])

    for (j, m) in enumerate(s)
        @test M[:, j] == m
    end
end

@testset "trees without branch lengths" begin
    s = splits(readnw("(T1,T2,(T3,T4));"))
    @test length(s) == 1
    @test isnan(only(collect(values(s.lengths))))
end
