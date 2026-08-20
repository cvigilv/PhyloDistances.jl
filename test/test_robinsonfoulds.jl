using Logging: Logging
using PhyloDistances
using PhyloDistances: isnormalized, normalizerinfo, splitinfo
using Random: Xoshiro
using Test

# Expected values produced by TreeDist 2.14.1 with ape 5.8.1:
#   RobinsonFoulds(t1, t2)  and  RobinsonFoulds(t1, t2, normalize = TRUE)
const RF_REFERENCE = [
    ("(A,B,(C,D));",             "(A,B,(C,D));",             0, 0.0),
    ("(A,B,(C,D));",             "(A,C,(B,D));",             2, 1.0),
    ("(A,B,((C,D),E));",         "(A,B,((C,E),D));",         2, 0.5),
    ("(A,B,(((C,D),E),F));",     "(A,B,(((C,E),D),F));",     2, 0.33333333333333331),
    ("(A,B,(((C,D),E),F));",     "(A,C,(((B,E),D),F));",     6, 1.0),
    ("((A,B),(C,D));",           "(A,B,(C,D));",             0, 0.0),
    ("(A,B,C,D,E);",             "(A,B,((C,D),E));",         2, 1.0),
    ("(A,B,(C,D),(E,F));",       "(A,B,(((C,D),E),F));",     3, 0.59999999999999998),
    ("(A,B,((C,D),(E,(F,G))));", "(A,B,((C,E),(D,(F,G))));", 4, 0.5),
]

@testset "RobinsonFoulds" begin
    @testset "matches the reference implementation" begin
        for (nw1, nw2, expected, expectednorm) in RF_REFERENCE
            t1, t2 = readnw(nw1), readnw(nw2)

            # A rooted input against an unrooted one warns; the value is what is at issue.
            got = @test_logs min_level = Logging.Error RobinsonFoulds()(t1, t2)
            @test got == expected

            gotnorm = @test_logs min_level = Logging.Error begin
                RobinsonFoulds(; normalize = true)(t1, t2)
            end
            @test gotnorm ≈ expectednorm
        end
    end

    @testset "two trees with no splits leave nothing to normalize by" begin
        # Both stars carry zero splits, so the reference reports 0/0.
        star = readnw("(A,B,C,D,E,F,G);")
        @test RobinsonFoulds()(star, star) == 0
        @test isnan(RobinsonFoulds(; normalize = true)(star, star))
    end

    @testset "identical trees are at distance zero" begin
        tree = randomtree(Xoshiro(1), 12)
        @test RobinsonFoulds()(tree, tree) == 0
        @test RobinsonFoulds(; normalize = true)(tree, tree) == 0
    end

    @testset "one rearrangement costs exactly two" begin
        rng = Xoshiro(4)
        for n in (5, 10, 25)
            tree = randomtree(rng, n)
            @test RobinsonFoulds()(tree, perturb(rng, tree, 1)) == 2
        end
    end

    @testset "binary trees sharing no split are 2(n-3) apart" begin
        # With n-3 splits each and none shared, the symmetric difference is everything.
        t1 = readnw("(A,B,(((C,D),E),F));")
        t2 = readnw("(A,C,(((B,E),D),F));")
        @test RobinsonFoulds()(t1, t2) == 2 * (6 - 3)
        @test RobinsonFoulds(; normalize = true)(t1, t2) == 1
    end

    @testset "branch lengths are ignored" begin
        @test RobinsonFoulds()(
            readnw("(A:1,B:1,(C:1,D:1):1);"), readnw("(A:9,B:8,(C:7,D:6):5);")
        ) == 0
    end

    @testset "symmetry" begin
        rng = Xoshiro(6)
        tree = randomtree(rng, 15)
        other = perturb(rng, tree, 8)
        @test RobinsonFoulds()(tree, other) == RobinsonFoulds()(other, tree)
    end

    @testset "normalization" begin
        t1, t2 = readnw("(A,B,((C,D),E));"), readnw("(A,B,((C,E),D));")

        @test normalizerinfo(RobinsonFoulds(), TreeDistConvention(), t1) == 2
        @test RobinsonFoulds(; normalize = max)(t1, t2) == 1.0
        @test RobinsonFoulds(; normalize = 4)(t1, t2) == 0.5
        @test !isnormalized(RobinsonFoulds())
    end

    @testset "result type" begin
        @test result_type(RobinsonFoulds(), Any, Any) === Int
        @test result_type(RobinsonFoulds(; normalize = true), Any, Any) === Float64

        trees = [randomtree(Xoshiro(i), 8) for i in 1:3]
        @test eltype(pairwise(RobinsonFoulds(), trees)) === Int
    end

    @testset "both conventions agree" begin
        t1, t2 = readnw("(A,B,((C,D),E));"), readnw("(A,B,((C,E),D));")
        @test RobinsonFoulds(; convention = :primary)(t1, t2) ==
              RobinsonFoulds(; convention = :treedist)(t1, t2)
    end

    @testset "mismatched taxa are rejected" begin
        @test_throws "trees span different taxa" begin
            RobinsonFoulds()(readnw("(A,B,(C,D));"), readnw("(A,B,(C,E));"))
        end
    end
end

@testset "WeightedRobinsonFoulds" begin
    @testset "identical topology and lengths give zero" begin
        tree = readnw("(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);")
        @test WeightedRobinsonFoulds()(tree, tree) == 0.0
    end

    @testset "shared splits contribute their length difference" begin
        # Same topology throughout; only C's branch differs, by 2.5.
        t1 = readnw("(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);")
        t2 = readnw("(A:1.0,B:2.0,(C:5.5,D:4.0):5.0);")
        @test WeightedRobinsonFoulds()(t1, t2) ≈ 2.5
    end

    @testset "a split in one tree only contributes its whole length" begin
        # {C,D} is worth 5 in the first tree and absent from the second, whose {B,C} is
        # worth 7. Pendant branches match, so the total is 5 + 7.
        t1 = readnw("(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);")
        t2 = readnw("(A:1.0,D:4.0,(B:2.0,C:3.0):7.0);")
        @test WeightedRobinsonFoulds()(t1, t2) ≈ 12.0
    end

    @testset "pendant branches count" begin
        # Identical topology, one pendant branch longer by 3.
        t1 = readnw("(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);")
        t2 = readnw("(A:4.0,B:2.0,(C:3.0,D:4.0):5.0);")
        @test WeightedRobinsonFoulds()(t1, t2) ≈ 3.0
    end

    @testset "a rooted tree matches its unrooted twin" begin
        # Unrooting sums the two root branches, which is inexact in floating point.
        rooted = readnw("((A:1.0,B:2.0):0.2,(C:3.0,D:4.0):0.4);")
        unrooted = readnw("(A:1.0,B:2.0,(C:3.0,D:4.0):0.6);")

        got = @test_logs min_level = Logging.Error WeightedRobinsonFoulds()(rooted, unrooted)
        @test got ≈ 0.0 atol = 1e-12
        @test got != 0.0
    end

    @testset "symmetry" begin
        t1 = readnw("(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);")
        t2 = readnw("(A:1.0,D:4.0,(B:2.0,C:3.0):7.0);")
        @test WeightedRobinsonFoulds()(t1, t2) == WeightedRobinsonFoulds()(t2, t1)
    end

    @testset "normalization is the combined branch length" begin
        t1 = readnw("(A:1.0,B:2.0,(C:3.0,D:4.0):5.0);")
        t2 = readnw("(A:1.0,D:4.0,(B:2.0,C:3.0):7.0);")

        @test normalizerinfo(WeightedRobinsonFoulds(), TreeDistConvention(), t1) ≈ 15.0
        @test normalizerinfo(WeightedRobinsonFoulds(), TreeDistConvention(), t2) ≈ 17.0
        @test WeightedRobinsonFoulds(; normalize = true)(t1, t2) ≈ 12.0 / 32.0
    end

    @testset "trees without branch lengths are rejected" begin
        @test_throws "every branch to have a length" begin
            WeightedRobinsonFoulds()(readnw("(A,B,(C,D));"), readnw("(A,B,(C,D));"))
        end
    end

    @testset "result type" begin
        @test result_type(WeightedRobinsonFoulds(), Any, Any) === Float64
    end
end

@testset "InfoRobinsonFoulds" begin
    @testset "hand-computed value for a five-taxon pair" begin
        # Each tree has exactly one non-trivial split ({C,D} vs {C,E}), and the two
        # differ, so nothing is shared and the distance is the sum of both trees' split
        # information.
        t1, t2 = readnw("(A,B,((C,D),E));"), readnw("(A,B,((C,E),D));")
        @test InfoRobinsonFoulds()(t1, t2) ≈ 2 * splitinfo(2, 5)
    end

    @testset "identical trees are at distance zero" begin
        tree = randomtree(Xoshiro(2), 12)
        @test InfoRobinsonFoulds()(tree, tree) == 0.0
        @test InfoRobinsonFoulds(; normalize = true)(tree, tree) == 0.0
    end

    @testset "a star tree costs the resolved tree's full split information" begin
        star = readnw("(A,B,C,D,E,F,G);")
        resolved = readnw("(A,B,(((C,D),E),F),G);")
        @test InfoRobinsonFoulds()(star, resolved) ≈
              normalizerinfo(InfoRobinsonFoulds(), TreeDistConvention(), resolved)
    end

    @testset "trees sharing no split sum their split information" begin
        t1 = readnw("(A,B,(((C,D),E),F));")
        t2 = readnw("(A,C,(((B,E),D),F));")
        expected = normalizerinfo(InfoRobinsonFoulds(), TreeDistConvention(), t1) +
                   normalizerinfo(InfoRobinsonFoulds(), TreeDistConvention(), t2)
        @test InfoRobinsonFoulds()(t1, t2) ≈ expected
        @test InfoRobinsonFoulds(; normalize = true)(t1, t2) == 1.0
    end

    @testset "a rooted tree matches its unrooted twin exactly" begin
        # Unlike WeightedRobinsonFoulds, no branch length is summed here, so the two
        # rootings agree bitwise rather than only to floating-point tolerance.
        rooted = readnw("((A,B),(C,D));")
        unrooted = readnw("(A,B,(C,D));")
        got = @test_logs min_level = Logging.Error InfoRobinsonFoulds()(rooted, unrooted)
        @test got == 0.0
    end

    @testset "branch lengths are ignored" begin
        @test InfoRobinsonFoulds()(
            readnw("(A:1,B:1,(C:1,D:1):1);"), readnw("(A:9,B:8,(C:7,D:6):5);")
        ) == 0.0
    end

    @testset "zero exactly where RobinsonFoulds is zero, positive where it isn't" begin
        rng = Xoshiro(9)
        for _ in 1:20
            t1 = randomtree(rng, rand(rng, 4:15))
            t2 = perturb(rng, t1, rand(rng, 0:3))
            rf = RobinsonFoulds()(t1, t2)
            irf = InfoRobinsonFoulds()(t1, t2)
            @test (rf == 0) == (irf == 0.0)
            @test irf >= 0.0
        end
    end

    @testset "symmetry" begin
        rng = Xoshiro(10)
        tree = randomtree(rng, 15)
        other = perturb(rng, tree, 8)
        @test InfoRobinsonFoulds()(tree, other) == InfoRobinsonFoulds()(other, tree)
    end

    @testset "normalization" begin
        t1, t2 = readnw("(A,B,((C,D),E));"), readnw("(A,B,((C,E),D));")
        expected = normalizerinfo(InfoRobinsonFoulds(), TreeDistConvention(), t1) +
                   normalizerinfo(InfoRobinsonFoulds(), TreeDistConvention(), t2)
        @test InfoRobinsonFoulds(; normalize = true)(t1, t2) ≈
              InfoRobinsonFoulds()(t1, t2) / expected
        @test !isnormalized(InfoRobinsonFoulds())
    end

    @testset "result type" begin
        @test result_type(InfoRobinsonFoulds(), Any, Any) === Float64
        @test result_type(InfoRobinsonFoulds(; normalize = true), Any, Any) === Float64
    end

    @testset "both conventions agree" begin
        t1, t2 = readnw("(A,B,((C,D),E));"), readnw("(A,B,((C,E),D));")
        @test InfoRobinsonFoulds(; convention = :primary)(t1, t2) ==
              InfoRobinsonFoulds(; convention = :treedist)(t1, t2)
    end

    @testset "mismatched taxa are rejected" begin
        @test_throws "trees span different taxa" begin
            InfoRobinsonFoulds()(readnw("(A,B,(C,D));"), readnw("(A,B,(C,E));"))
        end
    end
end
