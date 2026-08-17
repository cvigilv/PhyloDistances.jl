using PhyloDistances
using PhyloDistances: isrooted, nni!, taxonlabels
using Random: Xoshiro
using Test

"""Every internal node of a valid tree branches; only leaves are childless."""
function internalnodes(tree)
    return [
        n for n in PhyloDistances.AbstractTrees.PostOrderDFS(tree)
        if !isempty(PhyloDistances.AbstractTrees.children(n))
    ]
end

nchildren(node) = length(PhyloDistances.AbstractTrees.children(node))

@testset "randomtree" begin
    @testset "structure" begin
        for n in 3:12
            tree = randomtree(Xoshiro(n), n)

            @test length(taxonlabels(tree)) == n
            @test sort(taxonlabels(tree)) == sort(["T$i" for i in 1:n])

            # A degree-3 root marks the tree unrooted.
            @test !isrooted(tree)
            @test nchildren(tree) == 3

            # Binary everywhere below the root, so no node branches trivially.
            for node in internalnodes(tree)
                @test nchildren(node) >= 2
                node === tree || @test nchildren(node) == 2
            end
        end
    end

    @testset "binary, so it carries n-3 informative splits" begin
        for n in 4:12
            @test length(splits(randomtree(Xoshiro(n), n))) == n - 3
        end
    end

    @testset "branch lengths" begin
        tree = randomtree(Xoshiro(1), 8)
        s = splits(tree; trivial = true)
        @test all(!isnan, values(s.lengths))
        @test all(0 <= v < 1 for v in values(s.lengths))

        @testset "the generator is configurable" begin
            fixed = randomtree(Xoshiro(1), 6; branchlength = _ -> 2.5)
            lengths = values(splits(fixed; trivial = true).lengths)
            # Only the two branches either side of the root are summed into one split.
            @test all(v == 2.5 || v == 5.0 for v in lengths)
        end
    end

    @testset "labels are configurable" begin
        tree = randomtree(Xoshiro(1), 4; labels = ["w", "x", "y", "z"])
        @test sort(taxonlabels(tree)) == ["w", "x", "y", "z"]
    end

    @testset "reproducible under a fixed seed" begin
        a = randomtree(Xoshiro(42), 10)
        b = randomtree(Xoshiro(42), 10)
        index = taxonindex(a, b)

        @test collect(splits(a, index)) == collect(splits(b, index))
        @test splits(a, index).lengths == splits(b, index).lengths

        # Different seeds should generally give different topologies.
        c = randomtree(Xoshiro(7), 10)
        @test collect(splits(a, index)) != collect(splits(c, index))
    end

    @testset "rejects impossible requests" begin
        @test_throws "at least three taxa" randomtree(Xoshiro(1), 2)
        @test_throws "got 2 labels for 4 taxa" randomtree(Xoshiro(1), 4; labels = ["a", "b"])
    end
end

@testset "nni!" begin
    @testset "preserves the taxon set and tree validity" begin
        tree = randomtree(Xoshiro(3), 9)
        before = sort(taxonlabels(tree))

        nni!(Xoshiro(1), tree)

        @test sort(taxonlabels(tree)) == before
        @test nchildren(tree) == 3
        for node in internalnodes(tree)
            node === tree || @test nchildren(node) == 2
        end
        @test length(splits(tree)) == 9 - 3
    end

    @testset "moves exactly one split" begin
        # An interchange relocates a single subtree, so exactly one split is replaced and
        # the symmetric difference is always two, never more.
        rng = Xoshiro(17)
        for n in (5, 9, 20)
            tree = randomtree(rng, n)
            index = taxonindex(tree)
            base = splits(tree, index)

            for _ in 1:10
                @test length(symdiff(base, splits(perturb(rng, tree, 1), index))) == 2
            end
        end
    end

    @testset "changes the topology" begin
        # One interchange moves exactly one subtree, so at least one split must differ.
        tree = randomtree(Xoshiro(5), 8)
        index = taxonindex(tree)
        before = splits(tree, index)

        nni!(Xoshiro(2), tree)

        @test !isempty(symdiff(before, splits(tree, index)))
    end

    @testset "a tree with no internal branch has no interchange" begin
        @test_throws "no internal branch" nni!(Xoshiro(1), randomtree(Xoshiro(1), 3))
        @test_throws "at least four taxa" nni!(Xoshiro(1), readnw("(A:1,B:1,C:1);"))
    end
end

@testset "perturb" begin
    @testset "leaves the original alone" begin
        tree = randomtree(Xoshiro(11), 10)
        index = taxonindex(tree)
        before = collect(splits(tree, index))

        perturb(Xoshiro(1), tree, 5)

        @test collect(splits(tree, index)) == before
    end

    @testset "zero moves gives an equivalent tree" begin
        tree = randomtree(Xoshiro(11), 10)
        index = taxonindex(tree)
        @test isempty(symdiff(splits(tree, index), splits(perturb(Xoshiro(1), tree, 0), index)))
    end

    @testset "distance grows with the number of moves" begin
        # Individually a later move can undo an earlier one, so this holds on average
        # rather than for any single pair of trees.
        rng = Xoshiro(2024)
        tree = randomtree(rng, 30)
        index = taxonindex(tree)
        base = splits(tree, index)

        meandistance(k) = sum(
            length(symdiff(base, splits(perturb(rng, tree, k), index))) for _ in 1:20
        ) / 20

        distances = [meandistance(k) for k in (1, 5, 20)]
        @test issorted(distances)
        @test distances[1] < distances[end]
    end

    @testset "the taxon set survives many moves" begin
        tree = randomtree(Xoshiro(9), 12)
        moved = perturb(Xoshiro(4), tree, 50)

        @test sort(taxonlabels(moved)) == sort(taxonlabels(tree))
        @test length(splits(moved)) == 12 - 3
    end

    @testset "reproducible under a fixed seed" begin
        tree = randomtree(Xoshiro(11), 10)
        index = taxonindex(tree)
        a = perturb(Xoshiro(8), tree, 6)
        b = perturb(Xoshiro(8), tree, 6)

        @test collect(splits(a, index)) == collect(splits(b, index))
    end

    @test_throws "must be non-negative" perturb(Xoshiro(1), randomtree(Xoshiro(1), 5), -1)
end
