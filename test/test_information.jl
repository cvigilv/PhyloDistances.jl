using PhyloDistances
using PhyloDistances: clusteringentropy, jointentropy, log2rooted, log2unrooted,
    mutualinformation, SplitInfoTable, splitinfo
using Random
using Test

@testset "log2rooted" begin
    # (2m-3)!! for m = 0:7, by convention 1 for m <= 1.
    known = [1, 1, 1, 3, 15, 105, 945, 10395]
    for (m, expected) in enumerate(known)
        @test 2.0^log2rooted(m - 1) ≈ expected
    end
    @test_throws "negative" log2rooted(-1)
end

@testset "log2unrooted" begin
    # (2n-5)!! for n = 3:8, by convention 1 for n <= 2.
    known = [1, 1, 1, 1, 3, 15, 105, 945]
    for (n, expected) in enumerate(known)
        @test 2.0^log2unrooted(n - 1) ≈ expected
    end
    @test_throws "negative" log2unrooted(-1)

    # Finite and monotonically increasing well past any tip count this package targets.
    @test isfinite(log2unrooted(1000))
    @test log2unrooted(1000) > log2unrooted(999)
end

@testset "splitinfo" begin
    # A trivial split (one taxon against the rest) occurs in every tree on the same taxa,
    # so it carries no information, for any n.
    for n in (3, 5, 8, 100, 1000)
        @test splitinfo(1, n) ≈ 0.0 atol = 1e-9
        @test splitinfo(n - 1, n) ≈ 0.0 atol = 1e-9
    end

    # Hand-computed: splitinfo(2, 5) = log2unrooted(5) - log2rooted(2) - log2rooted(3)
    #               = log2(15) - log2(1) - log2(3) = log2(5).
    @test splitinfo(2, 5) ≈ log2(5)
    @test splitinfo(3, 5) ≈ log2(5)  # symmetric under k <-> n - k

    @test splitinfo(BitVector([0, 1, 1, 0, 0])) ≈ splitinfo(2, 5)

    @test_throws "0 <= k <= n" splitinfo(-1, 5)
    @test_throws "0 <= k <= n" splitinfo(6, 5)
end

@testset "SplitInfoTable" begin
    # The table is a cost optimization only, so every value must match the per-call form
    # exactly -- === rather than ≈, which also pins -0.0 against 0.0 and NaN against NaN.
    for n in (0, 1, 2, 3, 5, 8, 64, 200, 1000)
        table = SplitInfoTable(n)
        @test log2unrooted(table) === log2unrooted(n)
        @test all(k -> log2rooted(table, k) === log2rooted(k), 0:n)
        @test all(k -> splitinfo(table, k) === splitinfo(k, n), 0:n)
    end

    table = SplitInfoTable(5)
    @test splitinfo(table, BitVector([0, 1, 1, 0, 0])) === splitinfo(2, 5)
    @test string(table) == "SplitInfoTable(5)"

    @test_throws "negative" SplitInfoTable(-1)
    @test_throws "0 <= k <= n" splitinfo(table, -1)
    @test_throws "0 <= k <= n" splitinfo(table, 6)
    @test_throws "0 <= k <= n" log2rooted(table, 6)

    # A mask over the wrong taxon set would otherwise be answered for the wrong tree size.
    @test_throws DimensionMismatch splitinfo(table, BitVector([0, 1, 1, 0]))

    # Summing a real tree's splits is where the table is meant to be used, and where a
    # difference from the per-call form would show up as a drifting total rather than a
    # single wrong lookup.
    rng = Random.Xoshiro(20260820)
    for n in (6, 30, 200)
        s = splits(randomtree(rng, n))
        tbl = SplitInfoTable(n)
        @test sum(mask -> splitinfo(tbl, mask), s; init = 0.0) ===
            sum(splitinfo, s; init = 0.0)
    end
end

@testset "clusteringentropy" begin
    @test clusteringentropy(4, 8) ≈ 1.0  # an even split is exactly 1 bit
    @test clusteringentropy(2, 8) < 1.0  # an uneven split carries less

    # A taxon set with everything on one side is not a partition at all.
    @test clusteringentropy(0, 8) == 0.0
    @test clusteringentropy(8, 8) == 0.0
    @test clusteringentropy(0, 0) == 0.0

    @test clusteringentropy(3, 8) ≈ clusteringentropy(5, 8)  # symmetric under k <-> n - k
    @test clusteringentropy(BitVector([0, 1, 1, 0])) ≈ clusteringentropy(2, 4)

    @test_throws "0 <= k <= n" clusteringentropy(-1, 5)
    @test_throws "0 <= k <= n" clusteringentropy(6, 5)
end

@testset "mutualinformation" begin
    m1 = BitVector([0, 1, 1, 0, 0, 1, 0, 0])
    m2 = BitVector([0, 1, 0, 1, 0, 1, 0, 0])

    # A split's mutual information with itself is exactly its own entropy.
    @test mutualinformation(m1, m1) ≈ clusteringentropy(m1)

    # Symmetric in its two arguments, as any mutual information must be.
    @test mutualinformation(m1, m2) ≈ mutualinformation(m2, m1)

    # Bounded below by 0 and above by the smaller of the two marginal entropies.
    @test mutualinformation(m1, m2) >= 0.0
    @test mutualinformation(m1, m2) <= min(clusteringentropy(m1), clusteringentropy(m2))

    @test_throws DimensionMismatch mutualinformation(m1, BitVector([0, 1, 1]))

    # Nonnegativity is a property of the formula, not just these two hand-picked masks.
    rng = Random.Xoshiro(20260820)
    for _ in 1:2000
        n = rand(rng, 3:24)
        a, b = rand(rng, Bool, n), rand(rng, Bool, n)
        @test mutualinformation(a, b) >= -1e-9
    end
end

@testset "jointentropy" begin
    m1 = BitVector([0, 1, 1, 0, 0, 1, 0, 0])
    m2 = BitVector([0, 1, 0, 1, 0, 1, 0, 0])

    @test jointentropy(m1, m2) ≈
        clusteringentropy(m1) + clusteringentropy(m2) - mutualinformation(m1, m2)

    # Joint entropy is at least as large as either marginal entropy.
    @test jointentropy(m1, m2) >= clusteringentropy(m1) - 1e-9
    @test jointentropy(m1, m2) >= clusteringentropy(m2) - 1e-9

    @test jointentropy(m1, m1) ≈ clusteringentropy(m1)  # self-joint reduces to the entropy
end
