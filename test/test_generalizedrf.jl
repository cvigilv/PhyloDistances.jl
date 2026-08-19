using PhyloDistances
using PhyloDistances: _jaccardscore
using Random: Xoshiro
using Test

"""A caterpillar tree on the given taxon `order`, unrooted (root of degree 3)."""
function orderedcaterpillar(order)
    nw = "$(order[end]):1.0"
    for t in reverse(order[3:(end - 1)])
        nw = "($nw,$t:1.0):1.0"
    end
    return readnw("($(order[1]):1.0,$(order[2]):1.0,$nw);")
end

mk(ntaxa, idxs) = BitVector([i in idxs for i in 1:ntaxa])

@testset "_jaccardscore, hand-computed" begin
    ntaxa = 6

    @testset "identical splits score 1 regardless of k" begin
        a = mk(ntaxa, [3, 4])
        for k in (1, 2, 5)
            @test _jaccardscore(a, a, ntaxa; k, allowconflict = true) == 1.0
            @test _jaccardscore(a, a, ntaxa; k, allowconflict = false) == 1.0
        end
    end

    @testset "a nested (compatible but unequal) pair, worked by hand" begin
        # {T3,T4} <= {T3,T4,T5}: nesting makes the pair *compatible* — it does not make
        # the splits equal, and the score reflects that rather than jumping to 1.
        # a_and_b=2, a_and_B=0, A_and_b=1, A_and_B=3, giving Jaccard indices 2/3, 3/4, 0,
        # 1/6 for (a,b), (A,B), (a,B), (A,b) respectively; the score is
        # max(min(2/3,3/4), min(0,1/6)) = 2/3.
        a, b = mk(ntaxa, [3, 4]), mk(ntaxa, [3, 4, 5])
        @test _jaccardscore(a, b, ntaxa; k = 1, allowconflict = true) ≈ 2 / 3
        # Nesting is compatible (a ⊆ b), so forbidding conflict changes nothing here.
        @test _jaccardscore(a, b, ntaxa; k = 1, allowconflict = false) ≈ 2 / 3
        @test _jaccardscore(a, b, ntaxa; k = 1, allowconflict = true) ==
              _jaccardscore(b, a, ntaxa; k = 1, allowconflict = true)
    end

    @testset "a conflicting pair, worked by hand" begin
        # {T3,T4} vs {T4,T5}: neither is nested in the other or its complement.
        # |a∩b|=1, |a\b|=1, |b\a|=1, |neither|=3, giving Jaccard indices 1/3, 3/5, 1/5, 1/5
        # for (a,b), (A,B), (a,B), (A,b) respectively; the score is max(min(1/3,3/5),
        # min(1/5,1/5)) = 1/3.
        a, b = mk(ntaxa, [3, 4]), mk(ntaxa, [4, 5])
        @test _jaccardscore(a, b, ntaxa; k = 1, allowconflict = true) ≈ 1 / 3
        @test _jaccardscore(a, b, ntaxa; k = 1, allowconflict = false) == 0.0

        # Raising to k shrinks any score below 1 monotonically toward 0.
        @test _jaccardscore(a, b, ntaxa; k = 2, allowconflict = true) ≈ (1 / 3)^2
    end
end

@testset "identical trees maximize NyeSimilarity and vanish under JaccardRobinsonFoulds" begin
    tree = orderedcaterpillar(["T$i" for i in 1:8])
    n = length(splits(tree))

    @test NyeSimilarity()(tree, tree) == n
    @test JaccardRobinsonFoulds()(tree, tree) == 0.0
    @test JaccardRobinsonFoulds(; k = 3, allowconflict = false)(tree, tree) == 0.0
end

@testset "NyeSimilarity and JaccardRobinsonFoulds(k=1) are the same scoring, similarity vs distance" begin
    rng = Xoshiro(20260819)
    for _ in 1:20
        n = rand(rng, 6:15)
        t1 = randomtree(rng, n)
        t2 = perturb(rng, t1, rand(rng, 1:(2n)))
        index = taxonindex(t1, t2)
        n1, n2 = length(splits(t1, index)), length(splits(t2, index))

        nye = NyeSimilarity()(t1, t2)
        jrf = JaccardRobinsonFoulds(; k = 1, allowconflict = true)(t1, t2)
        @test jrf ≈ n1 + n2 - 2 * nye
    end
end

@testset "as k grows, JaccardRobinsonFoulds rises monotonically to RobinsonFoulds" begin
    # Only an identical pair of splits scores exactly 1 (the hand-computed cases above:
    # a nested-but-unequal pair falls short), so raising every other pair's score to a
    # growing power drives it to 0 and leaves only the exactly-shared splits contributing
    # — which is exactly what RobinsonFoulds counts. A higher exponent never scores a
    # mismatched pair *more* similar, so the distance cannot fall as k rises either.
    t1 = orderedcaterpillar(["T1", "T2", "T3", "T4", "T5", "T6", "T7"])
    t2 = orderedcaterpillar(["T1", "T2", "T4", "T3", "T6", "T5", "T7"])
    rf = RobinsonFoulds()(t1, t2)

    values = [JaccardRobinsonFoulds(; k)(t1, t2) for k in (1, 2, 5, 10, 30, 80, 200)]
    @test issorted(values)
    @test all(<=(rf), values)
    @test values[end] ≈ rf
end

@testset "allowconflict = false never scores lower than allowing conflicting matches" begin
    # A random pair (Xoshiro(555), n=7 from randomtree) whose optimal matching under
    # allowconflict = true actually uses a conflicting pair — a real exercise of the flag,
    # not a case where the two happen to coincide, which is common enough that most
    # hand-built examples do not exhibit it.
    t1 = readnw(
        "(T1:0.7630010824174914,(T2:0.9615006667557918,T7:0.1496603783858097)" *
        "0.0:0.4439412528792469,(T3:0.30625955035369534,(T4:0.6573482907349962," *
        "(T5:0.062030002357176595,T6:0.3656007970953302)0.0:0.6110904460076944)" *
        "0.0:0.03007992430872175)0.0:0.19215289428245164);"
    )
    t2 = readnw(
        "((T1:0.9331216533681328,T5:0.4090910016097077)0.0:0.29526737038348216," *
        "(T2:0.09291559672954408,((T4:0.4243357502882611,T7:0.4889904104608994)" *
        "0.0:0.25908297478850106,T6:0.283074948973815)0.0:0.8680507714687015)" *
        "0.0:0.2385858085469572,T3:0.817645617704037);"
    )

    allowed = JaccardRobinsonFoulds(; k = 1, allowconflict = true)(t1, t2)
    forbidden = JaccardRobinsonFoulds(; k = 1, allowconflict = false)(t1, t2)
    @test allowed ≈ 4.6 atol = 1e-9
    @test forbidden ≈ 6.2 atol = 1e-9
    @test allowed < forbidden
end

@testset "normalization" begin
    t1 = orderedcaterpillar(["T1", "T2", "T3", "T4", "T5", "T6", "T7"])
    t2 = orderedcaterpillar(["T1", "T2", "T4", "T3", "T6", "T5", "T7"])
    index = taxonindex(t1, t2)
    n1, n2 = length(splits(t1, index)), length(splits(t2, index))

    nye = NyeSimilarity()(t1, t2)
    @test NyeSimilarity(; normalize = true)(t1, t2) ≈ nye / ((n1 + n2) / 2)

    jrf = JaccardRobinsonFoulds()(t1, t2)
    @test JaccardRobinsonFoulds(; normalize = true)(t1, t2) ≈ jrf / (n1 + n2)
end

@testset "matched trees against a random collection are symmetric" begin
    rng = Xoshiro(1)
    trees = [randomtree(rng, 8) for _ in 1:4]
    for m in (NyeSimilarity(), JaccardRobinsonFoulds(; k = 2, allowconflict = false))
        for i in eachindex(trees), j in eachindex(trees)
            @test m(trees[i], trees[j]) ≈ m(trees[j], trees[i])
        end
    end
end
