# Vendored linear-assignment solver. Every generalized Robinson-Foulds metric reduces to
# an optimal matching between two trees' splits (see splitmatching.jl), and no maintained
# Julia package currently supplies the underlying assignment routine — hence vendoring it
# here rather than depending on one, per the project's design decisions.

"""
    _hungarian(cost::AbstractMatrix{<:Real}) -> (rowmatch, colmatch, total)

The minimum-cost assignment between the rows and columns of `cost`: a matching that pairs
every index on the smaller side to a distinct index on the larger side, minimizing the sum
of the matched entries.

`rowmatch[i]` is the column matched to row `i`, or `0` if row `i` goes unmatched (only
possible when there are more rows than columns); `colmatch` gives the same from the
column's side. `total` is the sum of `cost[i, rowmatch[i]]` over every matched row.

The classical shortest-augmenting-path assignment algorithm with vertex potentials (Kuhn
1955; Munkres 1957; the O(n²m) formulation), `O(min(n, m)² * max(n, m))` overall. Internal:
[`splitmatching`](@ref) is the documented entry point that turns a matching into a score.
"""
function _hungarian(cost::AbstractMatrix{<:Real})
    nr, nc = size(cost)
    if nr <= nc
        rowmatch, total = _hungarianwide(float.(cost))
        colmatch = zeros(Int, nc)
        for i in eachindex(rowmatch)
            rowmatch[i] == 0 && continue
            colmatch[rowmatch[i]] = i
        end
        return rowmatch, colmatch, total
    else
        colmatch, total = _hungarianwide(float.(permutedims(cost)))
        rowmatch = zeros(Int, nr)
        for j in eachindex(colmatch)
            colmatch[j] == 0 && continue
            rowmatch[colmatch[j]] = j
        end
        return rowmatch, colmatch, total
    end
end

# Requires n <= m; every row is matched to a distinct column. Column-indexed working
# arrays (`v`, `p`, `way`, `minv`, `used`) are sized `m + 1` and indexed `1:(m + 1)`,
# with index 1 standing in for the algorithm's usual sentinel column `0` — the classical
# presentation is 0-indexed throughout, and shifting by one is simpler than reproducing
# that indexing with OffsetArrays for a routine this small.
function _hungarianwide(cost::AbstractMatrix{Float64})
    n, m = size(cost)
    u = zeros(Float64, n)
    v = zeros(Float64, m + 1)
    p = zeros(Int, m + 1)
    way = zeros(Int, m + 1)

    for i in 1:n
        p[1] = i
        j0 = 1
        minv = fill(Inf, m + 1)
        used = falses(m + 1)
        while true
            used[j0] = true
            i0 = p[j0]
            delta = Inf
            j1 = 0
            for j in 2:(m + 1)
                used[j] && continue
                cur = cost[i0, j - 1] - u[i0] - v[j]
                if cur < minv[j]
                    minv[j] = cur
                    way[j] = j0
                end
                if minv[j] < delta
                    delta = minv[j]
                    j1 = j
                end
            end
            for j in 1:(m + 1)
                if used[j]
                    u[p[j]] += delta
                    v[j] -= delta
                else
                    minv[j] -= delta
                end
            end
            j0 = j1
            p[j0] != 0 || break
        end
        while true
            j1 = way[j0]
            p[j0] = p[j1]
            j0 = j1
            j0 != 1 || break
        end
    end

    rowmatch = zeros(Int, n)
    for j in 2:(m + 1)
        p[j] == 0 && continue
        rowmatch[p[j]] = j - 1
    end
    return rowmatch, -v[1]
end
