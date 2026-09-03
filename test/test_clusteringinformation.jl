using PhyloDistances
using PhyloDistances: _matchtotal, _mutualinformationscorematrix, clusteringentropy,
    mutualinformation, normalizerinfo
using Random: Xoshiro
using Test

"""A caterpillar tree with a prescribed taxon order and an unrooted degree-three root."""
function clusteringcaterpillar(order)
    nw = "$(order[end])"
    for taxon in reverse(order[3:(end - 1)])
        nw = "($nw,$taxon)"
    end
    return readnw("($(order[1]),$(order[2]),$nw);")
end

@testset "clustering information metrics" begin
    @testset "a one-split comparison reduces to the split-level quantities" begin
        t1 = readnw("(A,B,(C,D));")
        t2 = readnw("(A,C,(B,D));")
        index = taxonindex(t1, t2)
        s1, s2 = splits(t1, index), splits(t2, index)
        @test length(s1) == length(s2) == 1

        mi = mutualinformation(only(s1), only(s2))
        entropy = clusteringentropy(only(s1)) + clusteringentropy(only(s2))
        @test MutualClusteringInfo()(t1, t2) ≈ mi
        @test ClusteringInfoDistance()(t1, t2) ≈ entropy - 2 * mi
    end

    @testset "matches TreeDist 2.14.1 on a non-trivial matching" begin
        t1 = readnw("(A,B,((C,D),E));")
        t2 = readnw("(A,B,((C,E),D));")

        @test MutualClusteringInfo()(t1, t2) ≈ 0.99092368847664014 atol = 1.0e-12
        @test ClusteringInfoDistance()(t1, t2) ≈ 1.901955000865394 atol = 1.0e-12
        @test MutualClusteringInfo(; normalize = true)(t1, t2) ≈
            0.51028532972534479 atol = 1.0e-12
        @test ClusteringInfoDistance(; normalize = true)(t1, t2) ≈
            0.48971467027465521 atol = 1.0e-12
    end

    @testset "identical trees attain maximal similarity and zero distance" begin
        tree = randomtree(Xoshiro(14), 12)
        entropy = normalizerinfo(MutualClusteringInfo(), TreeDistConvention(), tree)

        @test MutualClusteringInfo()(tree, tree) == entropy
        @test MutualClusteringInfo(; normalize = true)(tree, tree) == 1.0
        @test ClusteringInfoDistance()(tree, tree) == 0.0
        @test ClusteringInfoDistance(; normalize = true)(tree, tree) == 0.0
    end

    @testset "unmatched splits contribute zero similarity and their full entropy to distance" begin
        star = readnw("(A,B,C,D,E,F,G);")
        resolved = readnw("(A,B,(((C,D),E),F),G);")
        entropy = normalizerinfo(ClusteringInfoDistance(), TreeDistConvention(), resolved)

        @test MutualClusteringInfo()(star, resolved) == 0.0
        @test ClusteringInfoDistance()(star, resolved) == entropy
        @test MutualClusteringInfo(; normalize = true)(star, resolved) == 0.0
        @test ClusteringInfoDistance(; normalize = true)(star, resolved) == 1.0
    end

    @testset "distance and similarity use the same optimal matching" begin
        rng = Xoshiro(20260821)
        for _ in 1:20
            n = rand(rng, 6:15)
            t1 = randomtree(rng, n)
            t2 = perturb(rng, t1, rand(rng, 1:n))
            h1 = normalizerinfo(MutualClusteringInfo(), TreeDistConvention(), t1)
            h2 = normalizerinfo(MutualClusteringInfo(), TreeDistConvention(), t2)
            mci = MutualClusteringInfo()(t1, t2)

            @test ClusteringInfoDistance()(t1, t2) ≈ h1 + h2 - 2 * mci
            @test MutualClusteringInfo(; normalize = true)(t1, t2) ≈
                mci / ((h1 + h2) / 2)
            @test ClusteringInfoDistance(; normalize = true)(t1, t2) ≈
                (h1 + h2 - 2 * mci) / (h1 + h2)
        end
    end

    @testset "exact split removal preserves the full assignment optimum" begin
        rng = Xoshiro(3614)
        sawshared = false
        for _ in 1:20
            n = rand(rng, 6:15)
            t1 = randomtree(rng, n)
            t2 = perturb(rng, t1, rand(rng, 1:n))
            index = taxonindex(t1, t2)
            s1, s2 = splits(t1, index), splits(t2, index)
            pairscore = _mutualinformationscorematrix(s1, s2, n)

            for j in axes(pairscore, 2), i in axes(pairscore, 1)
                @test pairscore[i, j] ≈ mutualinformation(s1.masks[i], s2.masks[j]) atol = 1.0e-14
            end

            fullscore = _matchtotal(pairscore; maximize = true).score
            reducedscore = MutualClusteringInfo()(t1, t2)
            @test reducedscore ≈ fullscore atol = 1.0e-12
            sawshared |= !isempty(intersect(s1, s2))
        end
        @test sawshared

        resolved = readnw("(A,B,(((C,D),E),F),G);")
        partial = readnw("(A,B,((C,D),E,F),G);")
        for (t1, t2) in ((resolved, partial), (partial, resolved))
            index = taxonindex(t1, t2)
            s1, s2 = splits(t1, index), splits(t2, index)
            @test length(intersect(s1, s2)) == min(length(s1), length(s2))
            pairscore = _mutualinformationscorematrix(s1, s2, length(index))
            fullscore = _matchtotal(pairscore; maximize = true).score
            @test MutualClusteringInfo()(t1, t2) ≈ fullscore atol = 1.0e-12
        end
    end

    @testset "symmetry and triangle inequality" begin
        rng = Xoshiro(1414)
        metric = ClusteringInfoDistance()
        for _ in 1:30
            n = rand(rng, 6:14)
            trees = [randomtree(rng, n) for _ in 1:3]
            d12 = metric(trees[1], trees[2])
            d21 = metric(trees[2], trees[1])
            d13 = metric(trees[1], trees[3])
            d23 = metric(trees[2], trees[3])

            @test d12 ≈ d21 atol = 1.0e-12
            @test d13 <= d12 + d23 + 1.0e-12
            @test d12 >= -1.0e-12
            @test MutualClusteringInfo()(trees[1], trees[2]) ≈
                MutualClusteringInfo()(trees[2], trees[1]) atol = 1.0e-12
        end
    end

    @testset "pairwise computes similarity diagonals" begin
        trees = [clusteringcaterpillar(circshift(["T$i" for i in 1:8], shift)) for shift in 0:2]
        similarities = pairwise(MutualClusteringInfo(), trees)
        distances = pairwise(ClusteringInfoDistance(), trees)

        @test similarities ≈ transpose(similarities)
        @test distances ≈ transpose(distances)
        @test all(i -> similarities[i, i] > 0, eachindex(trees))
        @test all(i -> distances[i, i] == 0, eachindex(trees))
    end

    @testset "configuration and input checks" begin
        t1 = readnw("(A,B,((C,D),E));")
        t2 = readnw("(A,B,((C,E),D));")

        for Comparison in (MutualClusteringInfo, ClusteringInfoDistance)
            primary = Comparison(; convention = :primary)(t1, t2)
            treedist = Comparison(; convention = :treedist)(t1, t2)
            @test primary == treedist
            @test result_type(Comparison(), Any, Any) === Float64
        end

        @test_throws "trees span different taxa" begin
            MutualClusteringInfo()(readnw("(A,B,(C,D));"), readnw("(A,B,(C,E));"))
        end
        @test_throws "trees span different taxa" begin
            ClusteringInfoDistance()(readnw("(A,B,(C,D));"), readnw("(A,B,(C,E));"))
        end
    end
end
