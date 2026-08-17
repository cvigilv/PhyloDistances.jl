using Distances: Distances
using OffsetArrays: OffsetVector
using PhyloDistances
using PhyloDistances: convention, isnormalized, issimilarity, normalization, normalizer,
    normalizerinfo, requiresrooted
using Test

# Applying a metric reconciles the rooting of its inputs, so these need real trees. The
# dummy metrics below reduce a tree to its leaf count, which keeps the expected value of
# every call obvious by inspection.

"""A star tree on `n` taxa: root of degree `n`, hence unrooted for `n >= 3`."""
star(n) = readnw("(" * join(("T$i:1.0" for i in 1:n), ",") * ");")

"""A tree on `n` taxa whose root has two children, hence rooted."""
rooted(n) = readnw("((" * join(("T$i:1.0" for i in 1:(n - 1)), ",") * "):1.0,T$n:1.0);")

nleaves(tree) = length(PhyloDistances.taxonlabels(tree))

"""Absolute difference; doubled under the primary convention, normalized by 10."""
struct Gap{C<:Convention,N} <: TreeMetric
    convention::C
    normalize::N
end
Gap(; convention = TreeDistConvention(), normalize = false) =
    Gap(Convention(convention), normalize)
PhyloDistances._compare(::Gap, ::TreeDistConvention, a, b) = abs(nleaves(a) - nleaves(b))
PhyloDistances._compare(::Gap, ::PrimaryConvention, a, b) = 2 * abs(nleaves(a) - nleaves(b))
PhyloDistances.normalizerinfo(::Gap, ::Convention, tree) = 5 * nleaves(tree)
Distances.result_type(m::Gap, ::Type, ::Type) = isnormalized(m) ? Float64 : Int

"""Implements only the TreeDist convention and defines no normalization."""
struct TreeDistOnly <: TreeMetric
    convention::Convention
    normalize::Any
end
TreeDistOnly(; convention = TreeDistConvention(), normalize = false) =
    TreeDistOnly(Convention(convention), normalize)
PhyloDistances._compare(::TreeDistOnly, ::TreeDistConvention, a, b) =
    abs(nleaves(a) - nleaves(b))

"""Largest on identical inputs, so a similarity rather than a distance."""
struct Closeness <: TreeSimilarity
    convention::Convention
    normalize::Any
end
Closeness(; convention = TreeDistConvention(), normalize = false) =
    Closeness(Convention(convention), normalize)
PhyloDistances._compare(::Closeness, ::TreeDistConvention, a, b) =
    100 - abs(nleaves(a) - nleaves(b))
Distances.result_type(::Closeness, ::Type, ::Type) = Int

"""Carries a metric-specific parameter and is defined on rooted trees."""
struct Scaled <: TreeMetric
    factor::Int
    convention::Convention
    normalize::Any
end
Scaled(factor; convention = TreeDistConvention(), normalize = false) =
    Scaled(factor, Convention(convention), normalize)
PhyloDistances._compare(m::Scaled, ::TreeDistConvention, a, b) =
    m.factor * abs(nleaves(a) - nleaves(b))
PhyloDistances.requiresrooted(::Scaled) = true

@testset "Convention" begin
    @test Convention(:treedist) === TreeDistConvention()
    @test Convention(:primary) === PrimaryConvention()

    # A Convention passes through, so constructors may accept either spelling.
    @test Convention(TreeDistConvention()) === TreeDistConvention()

    @test Symbol(TreeDistConvention()) === :treedist
    @test Symbol(PrimaryConvention()) === :primary

    @test_throws "unknown convention :nonesuch" Convention(:nonesuch)
    @test_throws "expected :treedist or :primary" Convention(:robinson)
end

@testset "type hierarchy" begin
    @test TreeMetric <: Distances.SemiMetric
    @test Gap() isa Distances.PreMetric

    # A similarity is deliberately outside the Distances hierarchy: SemiMetric's pairwise
    # takes the zero diagonal on faith, which would report zero self-similarity.
    @test !(TreeSimilarity <: Distances.PreMetric)
    @test Closeness() isa TreeComparison
    @test Gap() isa TreeComparison
end

@testset "traits" begin
    @test convention(Gap()) === TreeDistConvention()
    @test convention(Gap(; convention = :primary)) === PrimaryConvention()

    @test isnormalized(Gap()) === false
    @test isnormalized(Gap(; normalize = true)) === true

    @test requiresrooted(Gap()) === false
    @test requiresrooted(Scaled(3)) === true

    @test issimilarity(Gap()) === false
    @test issimilarity(Closeness()) === true
end

@testset "applying a metric" begin
    @testset "entry points agree" begin
        @test Gap()(star(3), star(10)) == 7
        @test evaluate(Gap(), star(3), star(10)) == 7
        @test evaluate(Closeness(), star(3), star(10)) == 93
    end

    @testset "symmetry" begin
        @test Gap()(star(3), star(10)) == Gap()(star(10), star(3))
    end

    @testset "convention is a field, not a keyword" begin
        @test Gap(; convention = :treedist)(star(3), star(10)) == 7
        @test Gap(; convention = :primary)(star(3), star(10)) == 14
        @test Gap(; convention = PrimaryConvention())(star(3), star(10)) == 14
        @test Gap()(star(3), star(10)) == Gap(; convention = :treedist)(star(3), star(10))
    end

    @testset "normalization is a field, not a keyword" begin
        # normalizerinfo is 5 * nleaves, and `true` sums it over both trees: 15 + 50 = 65.
        @test Gap(; normalize = true)(star(3), star(10)) == 7 / 65

        # Normalization applies to the value the convention produced.
        @test Gap(; convention = :primary, normalize = true)(star(3), star(10)) == 14 / 65
    end

    @testset "normalize accepts more than a flag" begin
        raw = Gap()(star(3), star(10))
        @test raw == 7

        @test normalization(Gap()) === false
        @test !isnormalized(Gap())
        @test isnormalized(Gap(; normalize = true))

        @testset "a number divides directly" begin
            @test Gap(; normalize = 10)(star(3), star(10)) == 0.7
            @test isnormalized(Gap(; normalize = 10))
        end

        @testset "a function combines the two trees' info" begin
            # normalizerinfo is 15 and 50 for these trees.
            @test Gap(; normalize = max)(star(3), star(10)) == 7 / 50
            @test Gap(; normalize = min)(star(3), star(10)) == 7 / 15
            @test Gap(; normalize = +)(star(3), star(10)) ==
                  Gap(; normalize = true)(star(3), star(10))
        end

        @testset "`true` is the metric's own scheme, not the number one" begin
            @test Gap(; normalize = true)(star(3), star(10)) != raw
            @test Gap(; normalize = 1)(star(3), star(10)) == raw
        end

        @test normalizerinfo(Gap(), TreeDistConvention(), star(3)) == 15
        @test normalizer(Gap(), TreeDistConvention(), star(3), star(10)) == 65
    end

    @testset "metric parameters live on the metric" begin
        @test Scaled(3)(rooted(3), rooted(7)) == 12
        @test Scaled(10)(rooted(3), rooted(7)) == 40
    end

    @testset "unimplemented convention reports the metric and the convention" begin
        @test TreeDistOnly()(star(3), star(10)) == 7
        @test_throws "TreeDistOnly does not implement the :primary convention" begin
            TreeDistOnly(; convention = :primary)(star(3), star(10))
        end
    end

    @testset "missing normalization is reported, not silently skipped" begin
        @test_throws "no normalization is defined for TreeDistOnly" begin
            TreeDistOnly(; normalize = true)(star(3), star(10))
        end
        @test_throws "no normalization is defined for TreeDistOnly" begin
            TreeDistOnly(; normalize = max)(star(3), star(10))
        end

        # A bare divisor needs nothing from the metric, so it still works.
        @test TreeDistOnly(; normalize = 2)(star(3), star(10)) == 3.5
    end
end

@testset "pairwise" begin
    trees = [star(3), star(6), star(13)]

    @testset "metric: agrees with repeated application" begin
        D = pairwise(Gap(), trees)
        for i in eachindex(trees), j in eachindex(trees)
            @test D[i, j] == Gap()(trees[i], trees[j])
        end
        @test D == [0 3 10; 3 0 7; 10 7 0]
    end

    @testset "metric: shape, symmetry and zero diagonal" begin
        D = pairwise(Gap(), trees)
        @test size(D) == (3, 3)
        @test D == D'
        @test all(iszero(D[i, i]) for i in axes(D, 1))
    end

    @testset "similarity: the diagonal is computed, not assumed" begin
        D = pairwise(Closeness(), trees)
        @test D == D'
        @test all(D[i, i] == 100 for i in axes(D, 1))
        @test D == [100 97 90; 97 100 93; 90 93 100]
    end

    @testset "two-collection form" begin
        @test pairwise(Closeness(), [star(3), star(6)], [star(13)]) == reshape([90, 93], 2, 1)
    end

    @testset "element type follows result_type" begin
        @test eltype(pairwise(Gap(), trees)) === Int
        @test eltype(pairwise(Gap(; normalize = true), trees)) === Float64
        @test eltype(pairwise(Closeness(), trees)) === Int
    end

    @testset "results are 1-based even for an offset input" begin
        # Distances.jl indexes observations by position rather than by key, and this
        # package follows it. Entries refer to positions in the collection.
        offset = OffsetVector(trees, -1:1)

        for m in (Gap(), Closeness())
            D = pairwise(m, offset)
            @test axes(D) == (Base.OneTo(3), Base.OneTo(3))
            @test D == pairwise(m, trees)
        end
    end
end
