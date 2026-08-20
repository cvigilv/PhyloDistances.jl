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

Jonker & Volgenant's (1987) shortest-augmenting-path algorithm: column reduction and
reduction transfer build a partial matching almost for free, augmenting row reduction
resolves most of what remains with only a row's two cheapest options, and a bucket-based
Dijkstra search settles whatever rows still need it. `O(n³)` worst case for a square
problem, `O(dim³)` here where `dim = max(n, m)` — internal:
[`splitmatching`](@ref) is the documented entry point that turns a matching into a score.
"""
function _hungarian(cost::AbstractMatrix{<:Real})
    nr, nc = size(cost)
    (nr == 0 || nc == 0) && return zeros(Int, nr), zeros(Int, nc), 0.0

    fcost = _asfloat(cost)
    dim = max(nr, nc)

    # A rectangular assignment reduces to a square one by padding the shorter side with
    # rows or columns of a single repeated value: however the `dim - min(nr, nc)` padding
    # entries are distributed among the padding rows/columns, they contribute the same
    # total, so the optimal square assignment always gives the real rows and columns their
    # best mutual matching and leaves the padding to soak up whatever is left over. This
    # holds for *any* single padding value, not just a large one — using the matrix's own
    # maximum keeps every entry on a comparable scale to the real costs, positive or
    # negative, rather than presuming a "large" constant is large enough.
    #
    # An already-square problem needs none of it, and `_jvlap` only reads its argument, so
    # the square case hands the caller's matrix straight through. Two trees on the same
    # taxa usually carry the same number of splits, which makes that the common case here
    # and the padded copy the largest single allocation the metric would otherwise make.
    padded = if nr == nc
        fcost
    else
        p = fill(maximum(fcost), dim, dim)
        @views p[1:nr, 1:nc] .= fcost
        p
    end

    # `_jvlap` reads its argument transposed, so it solves the problem whose rows are this
    # matrix's columns; its two results come back in the opposite order.
    colsol, rowsol = _jvlap(padded, dim)

    rowmatch = zeros(Int, nr)
    for i in 1:nr
        j = rowsol[i]
        j <= nc && (rowmatch[i] = j)
    end
    colmatch = zeros(Int, nc)
    for j in 1:nc
        i = colsol[j]
        i <= nr && (colmatch[j] = i)
    end

    total = 0.0
    for i in 1:nr
        rowmatch[i] == 0 && continue
        total += fcost[i, rowmatch[i]]
    end
    return rowmatch, colmatch, total
end

# `float.(cost)` always allocates a fresh array, even when `cost` is already `Float64` —
# wasteful for the common case of a cost matrix this package built itself. The algorithm
# only ever reads `cost`, so passing a `Float64` input through unchanged is safe.
_asfloat(cost::AbstractMatrix{Float64}) = cost
_asfloat(cost::AbstractMatrix{<:Real}) = float.(cost)

# The smallest value in column `j` of a `dim x dim` cost matrix and the (first, on ties)
# row that attains it. This is the one scan that runs against the storage order — see
# `_jvlap` on why the matrix is held transposed — and it is also the only one that visits
# every entry exactly once, rather than once per augmenting path.
function _jvcolmin(cost::AbstractMatrix{Float64}, j::Int, dim::Int)
    minval = cost[j, 1]
    imin = 1
    for i in 2:dim
        c = cost[j, i]
        if c < minval
            minval = c
            imin = i
        end
    end
    return minval, imin
end

# The two smallest values of `cost[j, i] - v[j]` over every column `j`, and the column
# attaining each. Augmenting row reduction (below) needs the runner-up as well as the
# best, to tell a column that is uniquely best for this row from one that is merely tied
# for best.
function _jvrowsubmin(cost::AbstractMatrix{Float64}, i::Int, v::Vector{Float64}, dim::Int)
    j1 = 1
    umin = cost[1, i] - v[1]
    usubmin = Inf
    j2 = 0
    for j in 2:dim
        h = cost[j, i] - v[j]
        if h < umin
            usubmin = umin
            j2 = j1
            umin = h
            j1 = j
        elseif h < usubmin
            usubmin = h
            j2 = j
        end
    end
    return umin, usubmin, j1, j2
end

# Whether `a` is less than `b` by more than floating-point noise, not literally less.
# Column potentials accumulate many small subtractions over the algorithm's run, so two
# columns that are genuinely tied for a row's best choice can end up differing by an ulp
# or two purely from rounding — and augmenting row reduction (below) uses "strictly less"
# to decide whether a column is *uniquely* best for a row. Without this margin, two rows
# tied at their best column can trade it back and forth forever: each swap "improves" the
# potential by an amount too small to actually change which column looks best next time,
# so the same false signal recurs. TreeDist's own solver guards the same comparison with
# `nontrivially_less_than`, which this mirrors — a genuine floating-point degeneracy, not
# merely a quantization artifact of its integer-scaled costs. The tolerance scales with
# the compared values' own magnitude so it stays meaningful whether costs are near zero
# or far from it.
_jvstrictlyless(a::Float64, b::Float64) = a < b - 8 * eps(max(abs(a), abs(b), 1.0))

"""
    _jvlap(cost::AbstractMatrix{Float64}, dim::Int) -> (rowsol, colsol)

The optimal assignment on a `dim x dim` cost matrix: `rowsol[i]` is the column matched to
row `i` and `colsol[j]` is the row matched to column `j`, a mutually consistent bijection
on `1:dim`. Every row and column is matched — there is no rectangular "leftover" here,
that is handled by `_hungarian`'s caller, which pads to reach this shape.

`cost` is held **transposed**: `cost[j, i]` is the cost of matching row `i` to column `j`.
Every phase after column reduction scans one row against all `dim` columns, so this puts
the scans down a column of storage, where Julia keeps consecutive elements. Reading a row
of a `Matrix` instead strides by `dim` elements, which at a thousand splits is a fresh
cache line per entry. An assignment problem and its transpose have the same solution with
the two sides exchanged, so a caller holding the natural orientation passes its matrix
unaltered and swaps the two results.

Jonker & Volgenant (1987), *A shortest augmenting path algorithm for dense and sparse
linear assignment problems*, Computing 38: 325–340. Column reduction greedily claims each
column's cheapest row, keeping only the tightest of several claims on one row; reduction
transfer tightens the columns that ended up claimed exactly once, using each such row's
second-best alternative; augmenting row reduction resolves most of what is still
unassigned using only a row's two cheapest columns, occasionally displacing a weaker
existing claim rather than searching further; and any row still unassigned after that is
settled by an explicit shortest-augmenting-path search (a Dijkstra variant that partitions
candidate columns into settled and frontier sets rather than rescanning all of them at
every step). `0` is this package's "unmatched" sentinel throughout, standing in for the
reference algorithm's `-1`.
"""
function _jvlap(cost::AbstractMatrix{Float64}, dim::Int)
    v = Vector{Float64}(undef, dim)
    rowsol = zeros(Int, dim)
    colsol = zeros(Int, dim)
    matches = zeros(Int, dim)

    # COLUMN REDUCTION. Processing columns from last to first is the reference algorithm's
    # own choice, not load-bearing here (any order visits every column exactly once and
    # reaches the same partial matching) — kept only to stay a faithful, checkable port.
    for j in dim:-1:1
        minval, imin = _jvcolmin(cost, j, dim)
        v[j] = minval
        matches[imin] += 1
        if matches[imin] == 1
            rowsol[imin] = j
            colsol[j] = imin
        elseif v[j] < v[rowsol[imin]]
            # A later (smaller-index) column claims this row more tightly than the
            # column that claimed it first; take the row over, freeing the loser.
            j1 = rowsol[imin]
            rowsol[imin] = j
            colsol[j] = imin
            colsol[j1] = 0
        else
            colsol[j] = 0
        end
    end

    # REDUCTION TRANSFER. A row claimed by exactly one column has no rival to arbitrate
    # against, so its column's potential can be tightened further using that row's
    # second-best alternative — free information the column-reduction pass above did not
    # use, since it only ever looked at each row's *best* option, never its second-best.
    freeunassigned = zeros(Int, dim)
    numfree = 0
    for i in 1:dim
        if matches[i] == 0
            numfree += 1
            freeunassigned[numfree] = i
        elseif matches[i] == 1
            j1 = rowsol[i]
            mincost = Inf
            for j in 1:dim
                j == j1 && continue
                rc = cost[j, i] - v[j]
                rc < mincost && (mincost = rc)
            end
            v[j1] -= mincost
        end
    end

    # AUGMENTING ROW REDUCTION, two passes (a third rarely helps and the reference
    # algorithm does not take one either). Each still-unassigned row looks at only its two
    # cheapest columns: if the best is *uniquely* best, the column's potential tightens
    # around it directly; otherwise, if the best column already belongs to another row,
    # this row takes it anyway and bumps the previous occupant back onto the free list —
    # a displacement chain that is still guaranteed to terminate in an improving
    # assignment, not a cycle, because each displacement strictly tightens some potential.
    loopcnt = 0
    while loopcnt < 2
        loopcnt += 1
        previous_numfree = numfree
        numfree = 0
        k = 1
        while k <= previous_numfree
            i = freeunassigned[k]
            k += 1
            umin, usubmin, j1, j2 = _jvrowsubmin(cost, i, v, dim)
            i0 = colsol[j1]
            strictly_less = _jvstrictlyless(umin, usubmin)
            if strictly_less
                v[j1] -= (usubmin - umin)
            elseif i0 != 0
                j1 = j2
                i0 = colsol[j2]
            end
            rowsol[i] = j1
            colsol[j1] = i
            if i0 != 0
                if strictly_less
                    # The displaced row might immediately find another uniquely-best
                    # column, so it is retried within *this* pass rather than deferred.
                    k -= 1
                    freeunassigned[k] = i0
                else
                    numfree += 1
                    freeunassigned[numfree] = i0
                end
            end
        end
    end

    # AUGMENT SOLUTION. Whatever rows augmenting row reduction could not place get a full
    # shortest-augmenting-path search each: `d` holds the current shortest known reduced
    # distance from `free_row` to each column, `pred` the row that distance was reached
    # through, and `cl[1:low-1]` / `cl[low:up-1]` / `cl[up:dim]` partition the columns into
    # settled, newly-settled-at-the-current-minimum, and frontier — so each round of
    # Dijkstra only rescans the frontier rather than every column.
    d = Vector{Float64}(undef, dim)
    pred = Vector{Int}(undef, dim)
    cl = Vector{Int}(undef, dim)

    for f in 1:numfree
        free_row = freeunassigned[f]
        for j in 1:dim
            d[j] = cost[j, free_row] - v[j]
            pred[j] = free_row
            cl[j] = j
        end

        unassignedfound = false
        endofpath = 0
        last = 0
        low = 1
        up = 1
        min_ = 0.0

        while !unassignedfound
            if up == low
                last = low - 1
                min_ = d[cl[up]]
                up += 1
                for k in up:dim
                    j = cl[k]
                    h = d[j]
                    if h <= min_
                        if h < min_
                            up = low
                            min_ = h
                        end
                        cl[k] = cl[up]
                        cl[up] = j
                        up += 1
                    end
                end
                for k in low:(up - 1)
                    if colsol[cl[k]] == 0
                        endofpath = cl[k]
                        unassignedfound = true
                        break
                    end
                end
            end

            if !unassignedfound
                j1 = cl[low]
                low += 1
                i = colsol[j1]
                h = cost[j1, i] - v[j1] - min_
                for k in up:dim
                    j = cl[k]
                    v2 = cost[j, i] - v[j] - h
                    if v2 < d[j]
                        pred[j] = i
                        if v2 == min_
                            if colsol[j] == 0
                                endofpath = j
                                unassignedfound = true
                                break
                            else
                                cl[k] = cl[up]
                                cl[up] = j
                                up += 1
                            end
                        end
                        d[j] = v2
                    end
                end
            end
        end

        for k in 1:last
            j1 = cl[k]
            v[j1] += d[j1] - min_
        end

        i = 0
        j1 = endofpath
        while i != free_row
            i = pred[j1]
            colsol[j1] = i
            j1, rowsol[i] = rowsol[i], j1
        end
    end

    return rowsol, colsol
end
