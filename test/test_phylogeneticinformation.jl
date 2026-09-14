using PhyloDistances
using PhyloDistances: _matchedphylogeneticinfo, _matchingsplitinfo, _matchtotal,
    _phylogeneticinfoscorematrix, _sharedphylogeneticinfo, normalizerinfo,
    SplitInfoTable, splitinfo
using Random: Xoshiro
using Test

phyloinfomask(n, members) = BitVector(i in members for i in 1:n)

function forcedphyloinfoexact(pairscore, s1, s2)
    unmatched1 = Int[]
    unmatched2 = collect(eachindex(s2.masks))
    fixedscore = 0.0
    for (i, mask) in pairs(s1.masks)
        j = findfirst(==(mask), s2.masks)
        if isnothing(j)
            push!(unmatched1, i)
        else
            fixedscore += pairscore[i, j]
            deleteat!(unmatched2, findfirst(==(j), unmatched2))
        end
    end
    (isempty(unmatched1) || isempty(unmatched2)) && return fixedscore
    return fixedscore + _matchtotal(
        copy(pairscore[unmatched1, unmatched2]); maximize = true
    ).score
end

@testset "phylogenetic information metrics" begin
    @testset "split-pair scores" begin
        table8 = SplitInfoTable(8)
        a = phyloinfomask(8, (3, 4, 5))
        nested = phyloinfomask(8, (3, 4, 5, 6))
        conflict = phyloinfomask(8, (4, 5, 6))

        @test _sharedphylogeneticinfo(table8, a, a) == splitinfo(table8, a)
        @test _sharedphylogeneticinfo(table8, a, .!a) == splitinfo(table8, a)
        # TreeDist::SplitSharedInformation(8, 3, 4), independently evaluated from
        # the number of binary trees consistent with both nested splits.
        @test _sharedphylogeneticinfo(table8, a, nested) ≈
            2.7224660244710979 atol = 1.0e-12
        @test _sharedphylogeneticinfo(table8, a, conflict) == 0.0
        @test _sharedphylogeneticinfo(table8, a, nested) ==
            _sharedphylogeneticinfo(table8, nested, a)

        # Smith's worked example: ABCDEF|GHIJ against ABCDEIJ|FGH. The most
        # informative consistent split is ABCDE|GH, a 5|2 split on seven taxa.
        m1 = phyloinfomask(10, 1:6)
        m2 = phyloinfomask(10, (1, 2, 3, 4, 5, 9, 10))
        @test _matchingsplitinfo(SplitInfoTable(10), m1, m2) ≈ splitinfo(5, 7)
        @test _matchingsplitinfo(SplitInfoTable(10), m1, m2) ==
            _matchingsplitinfo(SplitInfoTable(10), m2, m1)
    end

    @testset "exact pairs remain in the assignment problem" begin
        # An exact split need not belong to the best global assignment for either scorer.
        # These two fixtures make fixing every exact pair first strictly worse, guarding
        # against applying MutualClusteringInfo's exact-pair reduction here.
        spi1 = readnw(
            "(((((((T1,T11),T12),T8),(T9,T15)),T7),(T4,(T10,(T13,T14))))," *
            "(T2,T5),(T3,T6));"
        )
        spi2 = readnw(
            "((T3,T6),(T10,(T4,(T13,T14))),(T2,T5)," *
            "((T15,(T9,T7)),((T8,T1),(T12,T11))));"
        )
        msi1 = readnw(
            "(((T1,T12),((T5,T10),T8)),(T2,((T6,T15),T9))," *
            "((T3,(T14,T16)),(T4,((T7,T11),T13))));"
        )
        msi2 = readnw(
            "(((T1,T12),(T2,((T6,T15),T9))),((T5,T10),T8)," *
            "((T3,(T14,T16)),(T4,((T7,T11),T13))));"
        )

        for (scorer, t1, t2) in (
            (_sharedphylogeneticinfo, spi1, spi2),
            (_matchingsplitinfo, msi1, msi2),
        )
            index = taxonindex(t1, t2)
            s1, s2 = splits(t1, index), splits(t2, index)
            table = SplitInfoTable(length(index))
            pairscore = _phylogeneticinfoscorematrix(scorer, s1, s2, table)
            fullscore = _matchtotal(copy(pairscore); maximize = true).score
            @test fullscore > forcedphyloinfoexact(pairscore, s1, s2)
            @test _matchedphylogeneticinfo(scorer, s1, s2, table) ≈ fullscore
        end
    end

    @testset "incompatible one-split trees" begin
        t1 = readnw("(A,B,(C,D));")
        t2 = readnw("(A,C,(B,D));")
        information = 2 * splitinfo(2, 4)

        @test SharedPhylogeneticInfo()(t1, t2) == 0.0
        @test SharedPhylogeneticInfo(; normalize = true)(t1, t2) == 0.0
        @test MatchingSplitInfoDistance()(t1, t2) ≈ information
        @test MatchingSplitInfoDistance(; normalize = true)(t1, t2) ≈ 1.0
    end

    @testset "matches TreeDist 2.14.1 on a non-trivial matching" begin
        t1 = readnw("(A,B,((C,D),E));")
        t2 = readnw("(A,B,((C,E),D));")

        @test SharedPhylogeneticInfo()(t1, t2) ≈
            2.3219280948873613 atol = 1.0e-12
        @test SharedPhylogeneticInfo(; normalize = true)(t1, t2) ≈
            0.49999999999999972 atol = 1.0e-12
        @test MatchingSplitInfoDistance()(t1, t2) ≈
            2.9478623766648231 atol = 1.0e-12
        @test MatchingSplitInfoDistance(; normalize = true)(t1, t2) ≈
            0.31739380551401447 atol = 1.0e-12
    end

    @testset "self-comparisons and unmatched splits" begin
        tree = randomtree(Xoshiro(1515), 12)
        info = normalizerinfo(SharedPhylogeneticInfo(), TreeDistConvention(), tree)

        @test SharedPhylogeneticInfo()(tree, tree) ≈ info atol = 1.0e-12
        @test SharedPhylogeneticInfo(; normalize = true)(tree, tree) ≈ 1.0
        @test MatchingSplitInfoDistance()(tree, tree) ≈ 0.0 atol = 1.0e-12
        @test MatchingSplitInfoDistance(; normalize = true)(tree, tree) ≈ 0.0 atol = 1.0e-12

        star = readnw("(A,B,C,D,E,F,G);")
        resolved = readnw("(A,B,(((C,D),E),F),G);")
        resolvedinfo = normalizerinfo(
            MatchingSplitInfoDistance(), TreeDistConvention(), resolved
        )
        @test SharedPhylogeneticInfo()(star, resolved) == 0.0
        @test MatchingSplitInfoDistance()(star, resolved) == resolvedinfo
        @test SharedPhylogeneticInfo(; normalize = true)(star, resolved) == 0.0
        @test MatchingSplitInfoDistance(; normalize = true)(star, resolved) == 1.0

        @test SharedPhylogeneticInfo()(star, star) == 0.0
        @test MatchingSplitInfoDistance()(star, star) == 0.0
        @test isnan(SharedPhylogeneticInfo(; normalize = true)(star, star))
        @test isnan(MatchingSplitInfoDistance(; normalize = true)(star, star))
    end

    @testset "distance uses the same matching score" begin
        rng = Xoshiro(20260905)
        for _ in 1:20
            n = rand(rng, 6:15)
            t1 = randomtree(rng, n)
            t2 = perturb(rng, t1, rand(rng, 1:n))
            index = taxonindex(t1, t2)
            s1, s2 = splits(t1, index), splits(t2, index)
            table = SplitInfoTable(n)
            msi = _matchedphylogeneticinfo(_matchingsplitinfo, s1, s2, table)
            info1 = sum(mask -> splitinfo(table, mask), s1; init = 0.0)
            info2 = sum(mask -> splitinfo(table, mask), s2; init = 0.0)

            @test MatchingSplitInfoDistance()(t1, t2) ≈
                info1 + info2 - 2 * msi atol = 1.0e-12
            @test MatchingSplitInfoDistance(; normalize = true)(t1, t2) ≈
                (info1 + info2 - 2 * msi) / (info1 + info2) atol = 1.0e-12
        end
    end

    @testset "symmetry and pairwise behavior" begin
        rng = Xoshiro(1500)
        trees = [randomtree(rng, 10) for _ in 1:4]
        for comparison in (SharedPhylogeneticInfo(), MatchingSplitInfoDistance())
            for i in eachindex(trees), j in eachindex(trees)
                @test comparison(trees[i], trees[j]) ≈
                    comparison(trees[j], trees[i]) atol = 1.0e-12
            end
        end

        similarities = pairwise(SharedPhylogeneticInfo(), trees)
        distances = pairwise(MatchingSplitInfoDistance(), trees)
        @test similarities ≈ transpose(similarities)
        @test distances ≈ transpose(distances)
        @test all(i -> similarities[i, i] > 0, eachindex(trees))
        @test all(i -> distances[i, i] == 0, eachindex(trees))
    end

    @testset "configuration and input checks" begin
        t1 = readnw("(A,B,((C,D),E));")
        t2 = readnw("(A,B,((C,E),D));")

        for Comparison in (SharedPhylogeneticInfo, MatchingSplitInfoDistance)
            @test Comparison(; convention = :primary)(t1, t2) ==
                Comparison(; convention = :treedist)(t1, t2)
            @test result_type(Comparison(), Any, Any) === Float64
            @test_throws "trees span different taxa" Comparison()(
                readnw("(A,B,(C,D));"), readnw("(A,B,(C,E));")
            )
        end
    end
end
