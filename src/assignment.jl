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
        rowmatch, total = _hungarianwide(_asfloat(cost))
        colmatch = zeros(Int, nc)
        for i in eachindex(rowmatch)
            rowmatch[i] == 0 && continue
            colmatch[rowmatch[i]] = i
        end
        return rowmatch, colmatch, total
    else
        colmatch, total = _hungarianwide(_asfloat(permutedims(cost)))
        rowmatch = zeros(Int, nr)
        for j in eachindex(colmatch)
            colmatch[j] == 0 && continue
            rowmatch[colmatch[j]] = j
        end
        return rowmatch, colmatch, total
    end
end

# `float.(cost)` always allocates a fresh array, even when `cost` is already `Float64` —
# wasteful for the common case of a cost matrix this package built itself. The algorithm
# only ever reads `cost`, so passing a `Float64` input through unchanged is safe.
_asfloat(cost::AbstractMatrix{Float64}) = cost
_asfloat(cost::AbstractMatrix{<:Real}) = float.(cost)

# Requires n <= m; every row is matched to a distinct column. Column-indexed working
# arrays (`v`, `p`, `way`, `minv`, `used`) are sized `m + 1` and indexed `1:(m + 1)`,
# with index 1 standing in for the algorithm's usual sentinel column `0` — the classical
# presentation is 0-indexed throughout, and shifting by one is simpler than reproducing
# that indexing with OffsetArrays for a routine this small.
function _hungarianwide(cost::AbstractMatrix{Float64})
    n, m = size(cost)

    # Row reduction: a dual-feasible starting point — cost[i,j] - u[i] - v[j] >= 0 for
    # every i, j, since u[i] is literally row i's minimum and v starts at zero — cheap
    # (one O(n*m) pass) and a much better starting guess than potentials of all zero, since
    # it gives the search below a real chance of finding each row's own best column on the
    # first pass rather than discovering the same reduction lazily through many augmenting
    # steps. It stays valid for cost matrices with negative entries (this package always
    # negates similarities before minimizing) because duality requires only feasibility,
    # not nonnegativity of `cost` itself.
    #
    # Column reduction — the natural next step for the *square* assignment problem — is
    # deliberately not added on top of this. With n < m, only n of the m columns end up
    # matched, and reducing every column by its cross-row minimum can tie several columns
    # at reduced cost zero for one row purely because of what unrelated rows prefer, not
    # because those columns are equally good choices for *this* row once only the matched
    # columns are counted. The tie-break below (first zero found wins) can then lock in a
    # column that leaves a strictly better matching unreachable — caught by the brute-force
    # comparison in test/test_assignment.jl on a 3x6 case where it silently returned 5
    # instead of the true optimum 3.
    u = Vector{Float64}(undef, n)
    for i in 1:n
        rowmin = cost[i, 1]
        for j in 2:m
            c = cost[i, j]
            c < rowmin && (rowmin = c)
        end
        u[i] = rowmin
    end

    v = zeros(Float64, m + 1)

    # The order rows are processed in never changes which matching is optimal — each row's
    # search only ever augments the *current* partial matching by one more row, maintaining
    # dual feasibility and complementary slackness regardless of which row goes next — but
    # it changes how much work finding it costs. Rows with a low minimum tend to have a
    # narrow set of good columns and force long eviction chains if processed late, after
    # their best options are already claimed by less picky rows; processing the tightest
    # rows first avoids that. Measured on a ~1000-split matching, this roughly halved the
    # number of augmenting-path steps taken.
    order = sortperm(u)

    p = zeros(Int, m + 1)
    way = zeros(Int, m + 1)

    # Reused across every row rather than reallocated: the augmenting-path search below
    # runs once per row, and a matching's split count can reach the hundreds, so a fresh
    # allocation per row here was, before this fix, the single largest allocation source
    # for split-matching metrics on trees of any real size.
    #
    # `used` is a dense `Vector{Bool}`, not a `BitVector`: it is read and written many
    # times per row in the innermost loop, and `BitVector`'s bit-packed storage costs a
    # mask-and-shift on every access — measurable here since `used` is checked once per
    # candidate column, i.e. as often as the cost lookup itself. The memory this trades
    # away is a single `Bool` per column, negligible next to `cost`.
    minv = Vector{Float64}(undef, m + 1)
    used = Vector{Bool}(undef, m + 1)

    for i in order
        p[1] = i
        j0 = 1
        fill!(minv, Inf)
        fill!(used, false)
        while true
            used[j0] = true
            i0 = p[j0]
            delta = Inf
            j1 = 0
            # @inbounds on this loop and the one right after it covers every access in
            # both: `j` ranges over 2:(m + 1) or 1:(m + 1), matching v/minv/way/used's
            # declared size exactly, and `i0 = p[j0]` (here) and `p[j]` (in the next loop,
            # only read where `used[j]` is true) are always row indices in 1:n by the
            # algorithm's own invariant — p[1] is set to the current row above, and every
            # other value p ever holds is copied from another p[·] that was itself
            # established the same way (see the path-reconstruction loop below). Bounds
            # checking these accesses showed up as a real share of run time on a
            # ~1000-split matching in profiling, which is why this earns the annotation
            # rather than being added by reflex; the claim is exercised by the
            # brute-force and negative-cost correctness tests in test/test_assignment.jl,
            # which check thousands of matrices including rectangular ones where a wrong
            # `i0` or `p[j]` would show up as a wrong total, not necessarily a crash.
            @inbounds for j in 2:(m + 1)
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
            @inbounds for j in 1:(m + 1)
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

    # Read off the matched entries directly rather than from a potential-sum identity
    # (`-v[1]` in earlier versions of this function): every row is matched since n <= m,
    # and this way the total is correct regardless of how u and v were initialized, with
    # nothing downstream depending on a particular choice of starting potentials being
    # "the" canonical one.
    total = 0.0
    for i in 1:n
        total += cost[i, rowmatch[i]]
    end
    return rowmatch, total
end
