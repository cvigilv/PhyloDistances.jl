using Logging: Logging
using PhyloDistances
using PhyloDistances: isnormalized, normalizerinfo
using Random: Xoshiro
using Test

# The quartet distance read off the split sets instead of the path lengths: a tree resolves
# `ab|cd` exactly when one of its splits puts a and b on one side and c and d on the other.
# Slower and independent of the implementation, so a disagreement is informative.
function splitquartet(masks, a, b, c, d)
    for mask in masks
        mask[a] == mask[b] && mask[c] == mask[d] && mask[a] != mask[c] && return 1
        mask[a] == mask[c] && mask[b] == mask[d] && mask[a] != mask[b] && return 2
        mask[a] == mask[d] && mask[b] == mask[c] && mask[a] != mask[b] && return 3
    end
    return 0
end

function quartetdistance_bysplits(t1, t2)
    index = taxonindex(t1, t2)
    n = length(index)
    m1 = collect(splits(t1, index))
    m2 = collect(splits(t2, index))

    differing = 0
    for a in 1:(n - 3), b in (a + 1):(n - 2), c in (b + 1):(n - 1), d in (c + 1):n
        differing += splitquartet(m1, a, b, c, d) != splitquartet(m2, a, b, c, d)
    end
    return differing
end

@testset "QuartetDistance" begin
    @testset "leaf-to-leaf path lengths" begin
        tree = readnw("(A,B,(C,D));")
        # A and B hang off the root; C and D sit one edge further out.
        @test PhyloDistances._topologicalpaths(tree, taxonindex(tree)) == [
            0 2 3 3
            2 0 3 3
            3 3 0 2
            3 3 2 0
        ]
    end

    @testset "identical trees are at distance zero" begin
        tree = randomtree(Xoshiro(1), 12)
        @test QuartetDistance()(tree, tree) == 0
        @test QuartetDistance(; normalize = true)(tree, tree) == 0
    end

    @testset "four taxa admit one quartet" begin
        @test QuartetDistance()(readnw("(A,B,(C,D));"), readnw("(A,C,(B,D));")) == 1
        @test QuartetDistance()(readnw("(A,B,(C,D));"), readnw("(A,B,(C,D));")) == 0
    end

    @testset "one rearrangement of five taxa moves two quartets" begin
        # {A,C,D,E} and {B,C,D,E} are read as (A,E)(C,D) and (B,E)(C,D) by the first tree
        # and as (A,D)(C,E) and (B,D)(C,E) by the second; the other three agree.
        t1 = readnw("(A,B,((C,D),E));")
        t2 = readnw("(A,B,((C,E),D));")
        @test QuartetDistance()(t1, t2) == 2
    end

    @testset "a resolved tree is maximally far from the star tree" begin
        for (resolved, star, n) in (
            ("(A,B,(C,D));", "(A,B,C,D);", 4),
            ("(A,B,((C,D),E));", "(A,B,C,D,E);", 5),
            ("(A,B,(((C,D),E),F));", "(A,B,C,D,E,F);", 6),
        )
            @test QuartetDistance()(readnw(resolved), readnw(star)) == binomial(n, 4)
            @test QuartetDistance(; normalize = true)(readnw(resolved), readnw(star)) == 1
        end
    end

    @testset "a polytomy resolves only the quartets its splits separate" begin
        # The lone split {D,E} resolves exactly the three quartets holding both D and E.
        t1 = readnw("(A,B,C,(D,E));")
        t2 = readnw("(A,B,C,D,E);")
        @test QuartetDistance()(t1, t2) == 3
    end

    @testset "a rooted tree matches its unrooted twin" begin
        # The root is a node with one branch below it, which lengthens the paths across it
        # without changing which quartets the tree resolves.
        rooted = readnw("((A,B),(C,D));")
        unrooted = readnw("(A,B,(C,D));")
        got = @test_logs min_level = Logging.Error QuartetDistance()(rooted, unrooted)
        @test got == 0
    end

    @testset "branch lengths are ignored" begin
        @test QuartetDistance()(
            readnw("(A:1,B:1,(C:1,D:1):1);"), readnw("(A:9,B:8,(C:7,D:6):5);")
        ) == 0
    end

    @testset "agrees with a split-based enumeration" begin
        rng = Xoshiro(11)
        for n in (6, 7, 9), moves in (1, 3, 8)
            t1 = randomtree(rng, n)
            t2 = perturb(rng, t1, moves)
            @test QuartetDistance()(t1, t2) == quartetdistance_bysplits(t1, t2)
        end

        # Polytomies at several severities, which the random generator never produces.
        polytomous = readnw.([
            "(A,B,C,D,E,F,G);",
            "(A,B,C,(D,E,F,G));",
            "(A,B,(C,D),(E,F,G));",
            "(A,B,((C,D),E),(F,G));",
            "(A,B,(((C,D),E),F),G);",
        ])
        for t1 in polytomous, t2 in polytomous
            @test QuartetDistance()(t1, t2) == quartetdistance_bysplits(t1, t2)
        end
    end

    @testset "symmetry" begin
        rng = Xoshiro(6)
        t1 = randomtree(rng, 10)
        t2 = perturb(rng, t1, 5)
        @test QuartetDistance()(t1, t2) == QuartetDistance()(t2, t1)
    end

    @testset "normalization is the number of quartets" begin
        t1, t2 = readnw("(A,B,((C,D),E));"), readnw("(A,B,((C,E),D));")

        @test normalizerinfo(QuartetDistance(), TreeDistConvention(), t1) == 5
        @test QuartetDistance(; normalize = true)(t1, t2) == 2 / 5
        @test QuartetDistance(; normalize = max)(t1, t2) == 2 / 5
        @test QuartetDistance(; normalize = 4)(t1, t2) == 0.5
        @test !isnormalized(QuartetDistance())
    end

    @testset "fewer than four taxa leave nothing to compare" begin
        # No quartet exists, so the distance is zero and there is nothing to divide by.
        tree = readnw("(A,B,C);")
        @test QuartetDistance()(tree, tree) == 0
        @test isnan(QuartetDistance(; normalize = true)(tree, tree))
    end

    @testset "result type" begin
        @test result_type(QuartetDistance(), Any, Any) === Int
        @test result_type(QuartetDistance(; normalize = true), Any, Any) === Float64

        trees = [randomtree(Xoshiro(i), 8) for i in 1:3]
        @test eltype(pairwise(QuartetDistance(), trees)) === Int
    end

    @testset "both conventions agree" begin
        t1, t2 = readnw("(A,B,((C,D),E));"), readnw("(A,B,((C,E),D));")
        @test QuartetDistance(; convention = :primary)(t1, t2) ==
              QuartetDistance(; convention = :treedist)(t1, t2)
    end

    @testset "mismatched taxa are rejected" begin
        @test_throws "trees span different taxa" begin
            QuartetDistance()(readnw("(A,B,(C,D));"), readnw("(A,B,(C,E));"))
        end
    end
end
