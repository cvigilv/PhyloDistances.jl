using OffsetArrays: OffsetVector
using PhyloDistances
using PhyloDistances: issimilarity, normalizer, requiresrooted
using Test

# The interface performs no tree operations, so integers stand in for trees throughout:
# they make the expected value of every call obvious by inspection.

"""Absolute difference; doubled under the primary convention, normalized by 10."""
struct Gap <: TreeMetric end
PhyloDistances._compare(::Gap, ::TreeDistConvention, a, b) = abs(a - b)
PhyloDistances._compare(::Gap, ::PrimaryConvention, a, b) = 2 * abs(a - b)
PhyloDistances.normalizer(::Gap, ::Convention, a, b) = 10

"""Implements only the TreeDist convention and defines no normalization."""
struct TreeDistOnly <: TreeMetric end
PhyloDistances._compare(::TreeDistOnly, ::TreeDistConvention, a, b) = abs(a - b)

"""Carries its parameter as a field, and is a rooted similarity."""
struct Scaled <: TreeMetric
    factor::Int
end
PhyloDistances._compare(m::Scaled, ::TreeDistConvention, a, b) = m.factor * abs(a - b)
PhyloDistances.requiresrooted(::Scaled) = true
PhyloDistances.issimilarity(::Scaled) = true

@testset "Convention" begin
    @test Convention(:treedist) === TreeDistConvention()
    @test Convention(:primary) === PrimaryConvention()

    # A Convention passes through, so callers may use either spelling.
    @test Convention(TreeDistConvention()) === TreeDistConvention()

    @test Symbol(TreeDistConvention()) === :treedist
    @test Symbol(PrimaryConvention()) === :primary

    @test_throws "unknown convention :nonesuch" Convention(:nonesuch)
    @test_throws "expected :treedist or :primary" Convention(:robinson)
end

@testset "traits" begin
    @test requiresrooted(Gap()) === false
    @test issimilarity(Gap()) === false

    @test requiresrooted(Scaled(3)) === true
    @test issimilarity(Scaled(3)) === true
end

@testset "compare" begin
    @testset "symmetry" begin
        @test compare(Gap(), 3, 10) == compare(Gap(), 10, 3)
    end

    @testset "convention selection" begin
        @test compare(Gap(), 3, 10; convention = :treedist) == 7
        @test compare(Gap(), 3, 10; convention = :primary) == 14

        # Convention instances are accepted alongside the symbol shorthand.
        @test compare(Gap(), 3, 10; convention = PrimaryConvention()) == 14

        @test compare(Gap(), 3, 10) == compare(Gap(), 3, 10; convention = :treedist)
    end

    @testset "normalization" begin
        @test compare(Gap(), 3, 10; normalize = true) == 0.7

        # Normalization applies to the value the convention produced, not to the default.
        @test compare(Gap(), 3, 10; convention = :primary, normalize = true) == 1.4
    end

    @testset "metric parameters live on the metric" begin
        @test compare(Scaled(3), 1, 5) == 12
        @test compare(Scaled(10), 1, 5) == 40
    end

    @testset "unimplemented convention reports the metric and the convention" begin
        @test compare(TreeDistOnly(), 3, 10) == 7
        @test_throws "TreeDistOnly does not implement the :primary convention" begin
            compare(TreeDistOnly(), 3, 10; convention = :primary)
        end
    end

    @testset "missing normalization is reported, not silently skipped" begin
        @test_throws "no normalization is defined for TreeDistOnly" begin
            compare(TreeDistOnly(), 3, 10; normalize = true)
        end
    end

    @testset "metrics are applied through compare, not by calling them" begin
        @test !hasmethod(Gap(), Tuple{Int,Int})
    end
end

@testset "pairwise" begin
    trees = [0, 3, 10]

    @testset "agrees with repeated compare" begin
        D = pairwise(Gap(), trees)
        for i in eachindex(trees), j in eachindex(trees)
            @test D[i, j] == compare(Gap(), trees[i], trees[j])
        end
    end

    @testset "shape, symmetry and diagonal" begin
        D = pairwise(Gap(), trees)
        @test size(D) == (3, 3)
        @test D == D'
        @test all(iszero(D[i, i]) for i in axes(D, 1))
        @test D == [0 3 10; 3 0 7; 10 7 0]
    end

    @testset "keywords reach the metric" begin
        @test pairwise(Gap(), trees; convention = :primary) == 2 .* pairwise(Gap(), trees)
        @test pairwise(Gap(), trees; normalize = true) == pairwise(Gap(), trees) ./ 10
    end

    @testset "element type follows the metric" begin
        @test eltype(pairwise(Gap(), trees)) === Int
        @test eltype(pairwise(Gap(), trees; normalize = true)) === Float64
    end

    @testset "empty collection" begin
        D = pairwise(Gap(), Int[])
        @test size(D) == (0, 0)
        @test eltype(D) === Float64
    end

    @testset "axes of the result track the input" begin
        offset = OffsetVector(trees, -1:1)
        D = pairwise(Gap(), offset)

        @test axes(D) == (-1:1, -1:1)
        for i in axes(offset, 1), j in axes(offset, 1)
            @test D[i, j] == compare(Gap(), offset[i], offset[j])
        end

        # Same values as the 1-based call, just relabelled.
        @test parent(D) == pairwise(Gap(), trees)
    end
end
