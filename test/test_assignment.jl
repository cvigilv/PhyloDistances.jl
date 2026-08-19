using PhyloDistances: _hungarian
using Random
using Test

"""Call `f(perm)` for every injection from `1:k` into `1:n`, exhausting the assignment
problem's search space directly rather than through a permutation-generating dependency."""
function _eachinjection(f, n, k, chosen = Int[], used = falses(n))
    if length(chosen) == k
        f(chosen)
        return nothing
    end
    for i in 1:n
        used[i] && continue
        used[i] = true
        push!(chosen, i)
        _eachinjection(f, n, k, chosen, used)
        pop!(chosen)
        used[i] = false
    end
    return nothing
end

"""The minimum-cost assignment by exhaustive search, as the oracle `_hungarian` is checked
against. Matches every index on the smaller side to a distinct index on the larger side."""
function _bruteforceassignment(cost)
    nr, nc = size(cost)
    best = Ref(Inf)
    if nr <= nc
        _eachinjection(nc, nr) do perm
            best[] = min(best[], sum(cost[i, perm[i]] for i in 1:nr))
        end
    else
        _eachinjection(nr, nc) do perm
            best[] = min(best[], sum(cost[perm[j], j] for j in 1:nc))
        end
    end
    return best[]
end

@testset "empty and singleton matrices" begin
    @test _hungarian(zeros(0, 0)) == (Int[], Int[], -0.0)
    @test _hungarian(zeros(0, 3)) == (Int[], [0, 0, 0], -0.0)
    @test _hungarian(zeros(3, 0)) == ([0, 0, 0], Int[], -0.0)

    rowmatch, colmatch, total = _hungarian(reshape([5.0], 1, 1))
    @test rowmatch == [1] && colmatch == [1] && total == 5.0
end

@testset "hand-computed square case" begin
    # The diagonal (cost 1 + 1 = 2) beats the off-diagonal (2 + 2 = 4).
    rowmatch, colmatch, total = _hungarian([1.0 2.0; 2.0 1.0])
    @test rowmatch == [1, 2]
    @test colmatch == [1, 2]
    @test total == 2.0
end

@testset "a matching's reported total matches its cost" begin
    rng = Xoshiro(20260819)
    for _ in 1:500
        nr, nc = rand(rng, 1:8), rand(rng, 1:8)
        cost = Float64.(rand(rng, 1:6, nr, nc))
        rowmatch, colmatch, total = _hungarian(cost)

        matched = [(i, rowmatch[i]) for i in eachindex(rowmatch) if rowmatch[i] != 0]
        @test length(matched) == min(nr, nc)
        @test length(unique(last.(matched))) == length(matched)
        @test sum(cost[i, j] for (i, j) in matched; init = 0.0) ≈ total
        @test all(colmatch[j] == i for (i, j) in matched)
    end
end

@testset "agreement with brute-force search, including rectangular and tied costs" begin
    # A small integer range makes tied costs common, which is exactly where a matching
    # algorithm's tie-breaking is most likely to go wrong.
    rng = Xoshiro(8102026)
    for _ in 1:1000
        nr, nc = rand(rng, 1:7), rand(rng, 1:7)
        cost = Float64.(rand(rng, 1:5, nr, nc))
        _, _, total = _hungarian(cost)
        @test total ≈ _bruteforceassignment(cost)
    end

    # All-tied costs stress the case where every assignment is optimal.
    for n in 1:6, m in 1:6
        cost = fill(3.0, n, m)
        _, _, total = _hungarian(cost)
        @test total ≈ _bruteforceassignment(cost)
    end
end
