using PhyloDistances
using PhyloDistances: isrooted, requiresrooted, taxonlabel, taxonlabels
using Test

@testset "taxonlabel" begin
    tree = readnw("((A:1,B:2):3,C:4);")
    @test taxonlabel(first(PhyloDistances.AbstractTrees.Leaves(tree))) == "A"

    @test_throws "taxonlabel is not defined for Int64" taxonlabel(1)
end

@testset "taxonlabels" begin
    @test sort(taxonlabels(readnw("((A:1,B:2):3,C:4);"))) == ["A", "B", "C"]

    # A repeated label makes the correspondence between two trees ambiguous, so it is
    # rejected rather than silently resolved to whichever leaf is visited last.
    @test_throws "repeated taxon labels: \"A\"" taxonlabels(readnw("((A:1,A:2):3,C:4);"))
    @test_throws "must be unique within a tree" taxonlabels(readnw("((A:1,A:2):3,C:4);"))
    @test_throws "\"A\", \"B\"" taxonlabels(readnw("((A:1,A:2),(B:3,B:4),C:5);"))
end

@testset "TaxonIndex" begin
    index = taxonindex(readnw("((C:1,A:2):3,B:4);"))

    @testset "labels are sorted, so positions depend only on the taxon set" begin
        @test taxa(index) == ["A", "B", "C"]
        @test index["A"] == 1
        @test index["B"] == 2
        @test index["C"] == 3
    end

    @testset "leaf order in the Newick string does not matter" begin
        @test taxonindex(readnw("((C:1,A:2):3,B:4);")) ==
              taxonindex(readnw("((A:1,B:2):3,C:4);"))
    end

    @test length(index) == 3
    @test haskey(index, "A")
    @test !haskey(index, "Z")
    @test_throws KeyError index["Z"]

    @test occursin("3 taxa", sprint(show, index))
end

@testset "taxonindex over two trees" begin
    @testset "matching taxa" begin
        index = taxonindex(readnw("((A:1,B:2):3,C:4);"), readnw("((B:1,C:2):3,A:4);"))
        @test taxa(index) == ["A", "B", "C"]
    end

    @testset "mismatched taxa are named in both directions" begin
        t1 = readnw("((A:1,B:2):3,C:4);")
        t2 = readnw("((A:1,B:2):3,D:4);")

        @test_throws "trees span different taxa" taxonindex(t1, t2)
        @test_throws "only in the first tree: \"C\"" taxonindex(t1, t2)
        @test_throws "only in the second tree: \"D\"" taxonindex(t1, t2)
    end

    @testset "a strict subset names only the surplus side" begin
        t1 = readnw("((A:1,B:2):3,C:4);")
        t2 = readnw("(((A:1,B:2):3,C:4):5,D:6);")

        @test_throws "only in the first tree: none" taxonindex(t1, t2)
        @test_throws "only in the second tree: \"D\"" taxonindex(t1, t2)
    end
end

@testset "isrooted" begin
    # Newick carries no explicit rooting marker; the root's child count is the convention.
    @test isrooted(readnw("((A:1,B:2):3,C:4);"))
    @test isrooted(readnw("((A:1,B:2):3,(C:4,D:5):6);"))

    @test !isrooted(readnw("(A:1,B:2,C:4);"))
    @test !isrooted(readnw("(A:1,(B:2,C:3):4,D:5);"))

    @test_throws "cannot tell whether a tree is rooted" isrooted(readnw("(A:1);"))
    @test_throws "root with 1 child" isrooted(readnw("(A:1);"))
    @test_throws "at least three taxa" isrooted(readnw("(A:1);"))
end

# Metrics used only to drive the rooting reconciliation; the value they return is
# irrelevant here.
struct NeedsUnrooted <: TreeMetric
    convention::Convention
    normalize::Any
end
NeedsUnrooted() = NeedsUnrooted(TreeDistConvention(), false)
PhyloDistances._compare(::NeedsUnrooted, ::TreeDistConvention, a, b) = 0.0

struct NeedsRooted <: TreeMetric
    convention::Convention
    normalize::Any
end
NeedsRooted() = NeedsRooted(TreeDistConvention(), false)
PhyloDistances._compare(::NeedsRooted, ::TreeDistConvention, a, b) = 0.0
PhyloDistances.requiresrooted(::NeedsRooted) = true

@testset "rooting reconciliation" begin
    rooted = readnw("((A:1,B:2):3,C:4);")
    unrooted = readnw("(A:1,B:2,C:4);")

    @testset "matching rooting is silent" begin
        @test_logs NeedsUnrooted()(unrooted, unrooted)
        @test_logs NeedsRooted()(rooted, rooted)
    end

    @testset "a rooted tree given to an unrooted metric warns and proceeds" begin
        # Defensible without asking the user: the root's two child branches induce the
        # same bipartition, so discarding the root leaves a well-defined unrooted tree.
        @test_logs (:warn, r"NeedsUnrooted is defined on unrooted trees") begin
            NeedsUnrooted()(rooted, unrooted)
        end
        @test_logs (:warn, r"the first tree is rooted") NeedsUnrooted()(rooted, unrooted)
        @test_logs (:warn, r"the second tree is rooted") NeedsUnrooted()(unrooted, rooted)

        # The warning states the conversion applied, which is what makes warning rather
        # than throwing acceptable here.
        @test_logs (:warn, r"root position is ignored") NeedsUnrooted()(rooted, unrooted)

        @test (@test_logs (:warn,) NeedsUnrooted()(rooted, unrooted)) == 0.0
    end

    @testset "an unrooted tree given to a rooted metric throws" begin
        # No canonical rooting exists, so the choice belongs to the caller.
        @test_throws "NeedsRooted is defined on rooted trees" NeedsRooted()(unrooted, rooted)
        @test_throws "the first tree is unrooted" NeedsRooted()(unrooted, rooted)
        @test_throws "the second tree is unrooted" NeedsRooted()(rooted, unrooted)
        @test_throws "root it explicitly before comparing" NeedsRooted()(unrooted, rooted)
    end
end
