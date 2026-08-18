using Logging: Logging
using PhyloDistances
using Test

# The fixture's format is defined once, and `validation/fixture.jl` reads it through the
# same code to regenerate the file from TreeDist and Quartet.
include(joinpath(@__DIR__, "fixtures", "read.jl"))

const FIXTURE = readfixture()

# `NaN` agrees only with `NaN`, which is what a normalized distance gives when neither
# tree carries a split and the divisor is zero. Everything else is compared bitwise: the
# committed values are those of the reference implementations, reached by the same IEEE
# operations, so a tolerance would only hide a real difference.
agrees(got::Integer, expected::Integer) = got == expected
agrees(got::AbstractFloat, expected::AbstractFloat) =
    isnan(expected) ? isnan(got) : got === expected

@testset "reference values" begin
    @test length(FIXTURE) == 9
    @test allunique(row.case for row in FIXTURE)

    # One case is a rooted tree compared under an unrooted metric, which warns by design;
    # the values are what is at issue here.
    for row in FIXTURE
        @testset "$(row.case)" begin
            t1, t2 = readnw(row.newick1), readnw(row.newick2)

            @test agrees(@test_logs(
                min_level = Logging.Error, RobinsonFoulds()(t1, t2)
            ), row.rf)
            @test agrees(@test_logs(
                min_level = Logging.Error, RobinsonFoulds(; normalize = true)(t1, t2)
            ), row.rf_normalized)

            @test agrees(@test_logs(
                min_level = Logging.Error, QuartetDistance()(t1, t2)
            ), row.quartet)
            @test agrees(@test_logs(
                min_level = Logging.Error, QuartetDistance(; normalize = true)(t1, t2)
            ), row.quartet_normalized)

            # Both metrics are symmetric, and neither depends on argument order for its
            # rooting reconciliation.
            @test agrees(@test_logs(
                min_level = Logging.Error, RobinsonFoulds()(t2, t1)
            ), row.rf)
            @test agrees(@test_logs(
                min_level = Logging.Error, QuartetDistance()(t2, t1)
            ), row.quartet)
        end
    end

    @testset "the rooted case warns and says so" begin
        row = only(r for r in FIXTURE if r.case == "rooted-vs-unrooted-twin")
        t1, t2 = readnw(row.newick1), readnw(row.newick2)
        @test PhyloDistances.isrooted(t1)
        @test !PhyloDistances.isrooted(t2)
        @test (@test_logs (:warn, r"rooted") RobinsonFoulds()(t1, t2)) == row.rf
        @test (@test_logs (:warn, r"rooted") QuartetDistance()(t1, t2)) == row.quartet
    end

    @testset "a malformed file is rejected rather than misread" begin
        path = tempname()
        try
            write(path, "case\tnewick1\n" * "a\t(A,B,(C,D));\n")
            @test_throws "expected" readfixture(path)

            write(path, join(FIXTURE_COLUMNS, '\t') * "\n" * "short\t(A,B,(C,D));\n")
            @test_throws "expected 8" readfixture(path)

            write(path, "# nothing but a comment\n")
            @test_throws "no header line" readfixture(path)
        finally
            rm(path; force = true)
        end
    end
end

@testset "all-pairs over a tree collection" begin
    # Four resolvings of six taxa and the star that resolves none of them. Trees 1, 2 and
    # 4 are binary and carry three splits; tree 3 has a polytomy at the root and carries
    # two, resolving 11 of the 15 quartets.
    trees = readnw.([
        "(A,B,((C,D),(E,F)));",
        "(A,B,(((C,D),E),F));",
        "(A,B,(C,D),(E,F));",
        "(A,C,(((B,E),D),F));",
        "(A,B,C,D,E,F);",
    ])
    n = length(trees)

    @testset "$(nameof(typeof(metric)))" for metric in
                                             (RobinsonFoulds(), QuartetDistance())
        D = pairwise(metric, trees)

        @test size(D) == (n, n)
        @test eltype(D) === Int
        @test all(iszero, D[i, i] for i in 1:n)
        @test D == transpose(D)
        for j in 1:n, i in 1:n
            @test D[i, j] == metric(trees[i], trees[j])
        end
    end

    @testset "the star is the far corner" begin
        rf = pairwise(RobinsonFoulds(), trees)
        quartet = pairwise(QuartetDistance(), trees)

        # Sharing no split with the star, each tree is its own split count away.
        @test rf[5, :] == [3, 3, 2, 3, 0]
        # The star resolves nothing, so every quartet either tree resolves is a difference.
        @test quartet[5, :] == [15, 15, 11, 15, 0]
        @test quartet[5, 1] == binomial(6, 4)
    end

    @testset "Robinson-Foulds obeys the triangle inequality" begin
        D = pairwise(RobinsonFoulds(), trees)
        for k in 1:n, j in 1:n, i in 1:n
            @test D[i, j] <= D[i, k] + D[k, j]
        end
    end

    @testset "normalizing leaves a matrix of fractions" begin
        for metric in (RobinsonFoulds(; normalize = true),
                       QuartetDistance(; normalize = true))
            D = pairwise(metric, trees)
            @test eltype(D) === Float64
            @test all(iszero, D[i, i] for i in 1:n)
            @test all(x -> 0 <= x <= 1, D)
        end
    end
end
