# Phylogenetic-information comparisons. Both metrics maximize a one-to-one matching of
# non-trivial splits, but they assign different information scores to each candidate pair.

# log2 of the number of unrooted binary trees on `k` leaves, read from a table built for a
# larger taxon set. SplitInfoTable stores the rooted recurrence at every size, so this stays
# constant-time without allocating one table per candidate split.
function _log2unrooted(table::SplitInfoTable, k::Integer)
    n = _ntaxa(table)
    (0 <= k <= n) || throw(ArgumentError(
        "log2unrooted needs 0 <= k <= n, got k=$k from a table built for n=$n"
    ))
    k <= 2 && return 0.0
    return table.rooted[k + 1] - log2(2k - 3)
end

# The logarithm of the number of trees consistent with either orientation of two nested
# split sides, after the terms for each split separately have been cancelled. This is
# TreeDist's `one_overlap` quantity.
function _oneoverlap(table::SplitInfoTable, a::Integer, b::Integer)
    n = _ntaxa(table)
    (0 <= a <= n && 0 <= b <= n) || throw(ArgumentError(
        "overlap sizes must lie between 0 and $n, got $a and $b"
    ))
    a == b && return table.rooted[a + 1] + table.rooted[n - a + 1]

    lo, hi = minmax(a, b)
    return table.rooted[hi + 1] + table.rooted[n - lo + 1] -
           table.rooted[hi - lo + 2]
end

# Shared phylogenetic information from the sufficient counts of two splits. Compatible
# splits have at least one empty region in their four-way Venn diagram. Conflicting splits
# have all four regions populated and share no phylogenetic information under this measure.
function _sharedphylogeneticinfo(
    table::SplitInfoTable, aandb::Integer, na::Integer, nb::Integer
)
    n = _ntaxa(table)
    nA = n - na
    aandB = na - aandb
    Aandb = nb - aandb
    AandB = nA - Aandb
    minimum((aandb, aandB, Aandb, AandB)) >= 0 || throw(ArgumentError(
        "invalid split contingency counts: intersection=$aandb, sizes=$na and $nb, n=$n"
    ))

    # Equal and complementary masks describe the same split. Use splitinfo's arithmetic
    # path so self-similarity and the per-tree normalizer agree exactly.
    ((aandB == 0 && Aandb == 0) || (aandb == 0 && AandB == 0)) &&
        return splitinfo(table, na)

    overlap = if aandb == 0
        _oneoverlap(table, na, n - nb)
    elseif aandb == nb || aandb == na
        _oneoverlap(table, na, nb)
    elseif na + nb - aandb == n
        _oneoverlap(table, na, n - nb)
    else
        return 0.0
    end
    return log2unrooted(table) - overlap
end

function _sharedphylogeneticinfo(
    table::SplitInfoTable, a::BitVector, b::BitVector
)
    n = _ntaxa(table)
    (length(a) == n && length(b) == n) || throw(DimensionMismatch(
        "shared phylogenetic information needs two masks over the table's $n taxa"
    ))
    return _sharedphylogeneticinfo(table, _countand(a, b), count(a), count(b))
end

# Information in a split of an `x`-taxon subset. A side of size zero or one is present in
# every tree on that subset and therefore carries no information.
function _subsetsplitinfo(table::SplitInfoTable, y::Integer, x::Integer)
    (0 <= y <= x <= _ntaxa(table)) || throw(ArgumentError(
        "subset split information needs 0 <= y <= x <= $(_ntaxa(table)), got y=$y, x=$x"
    ))
    (y <= 1 || x - y <= 1) && return 0.0
    return _log2unrooted(table, x) - table.rooted[y + 1] - table.rooted[x - y + 1]
end

# Matching split information from the same four contingency counts. The two candidate
# splits separate taxa that agree between the input splits, or taxa that differ. The more
# informative candidate supplies the pair score.
function _matchingsplitinfo(
    table::SplitInfoTable, aandb::Integer, na::Integer, nb::Integer
)
    n = _ntaxa(table)
    ndifferent = na + nb - 2aandb
    nsame = n - ndifferent
    naonly = na - aandb

    if nsame == 0 || nsame == aandb
        return _subsetsplitinfo(table, naonly, ndifferent)
    elseif ndifferent == 0 || ndifferent == naonly
        return _subsetsplitinfo(table, aandb, nsame)
    end

    return max(
        _subsetsplitinfo(table, aandb, nsame),
        _subsetsplitinfo(table, naonly, ndifferent),
    )
end

function _matchingsplitinfo(table::SplitInfoTable, a::BitVector, b::BitVector)
    n = _ntaxa(table)
    (length(a) == n && length(b) == n) || throw(DimensionMismatch(
        "matching split information needs two masks over the table's $n taxa"
    ))
    return _matchingsplitinfo(table, _countand(a, b), count(a), count(b))
end

# Build either information score matrix from packed masks. Marginal split sizes belong to
# one row or column, leaving only the intersection population to count per matrix cell.
function _phylogeneticinfoscorematrix(
    scorer, s1::Splits, s2::Splits, table::SplitInfoTable
)
    _checksharedindex(s1, s2)
    _ntaxa(table) == length(taxonindex(s1)) || throw(DimensionMismatch(
        "the information table and split sets must span the same number of taxa"
    ))

    w1, w2 = _packmasks(s1.masks), _packmasks(s2.masks)
    na = [count(mask) for mask in s1.masks]
    nb = [count(mask) for mask in s2.masks]
    pairscore = Matrix{Float64}(undef, length(s1), length(s2))
    nwords = size(w1, 1)

    for j in axes(pairscore, 2)
        nbj = nb[j]
        for i in axes(pairscore, 1)
            pairscore[i, j] = scorer(
                table, _countand(w1, w2, i, j, nwords), na[i], nbj
            )
        end
    end
    return pairscore
end

function _matchedphylogeneticinfo(scorer, s1::Splits, s2::Splits, table::SplitInfoTable)
    (isempty(s1) || isempty(s2)) && return 0.0
    pairscore = _phylogeneticinfoscorematrix(scorer, s1, s2, table)
    return _matchtotal(pairscore; maximize = true).score
end

_treephylogeneticinfo(tree) = begin
    index = taxonindex(tree)
    _splitwiseinfo(SplitInfoTable(length(index)), splits(tree, index))
end

"""
    SharedPhylogeneticInfo(; convention = :treedist, normalize = false)

The shared phylogenetic information similarity of Smith (2020), in bits. Each pair of
compatible splits is scored by the information they share: the negative base-two logarithm
of the fraction of unrooted binary trees consistent with either split, corrected for the
information in each split separately. Incompatible splits score exactly zero. The tree
similarity is the largest total over a one-to-one matching of non-trivial splits, with
excess splits left unmatched.

The self-similarity of a tree is its total splitwise information, [`splitinfo`](@ref)
summed over its non-trivial splits. With `normalize = true`, the result is divided by the
mean splitwise information of the two trees, following TreeDist. A function such as `max`
or `min` combines the two per-tree totals instead, and a number is used directly as the
divisor.

Both conventions currently compute the same definition. Computing all split-pair scores
takes `O(s₁ s₂ ceil(n / 64))` time, followed by a dense linear assignment solve, where `n`
is the taxon count and `s₁`, `s₂` are the split counts.

Smith, M.R. (2020). *Information theoretic Generalized Robinson-Foulds metrics for
comparing phylogenetic trees.* Bioinformatics 36(20): 5007–5013.
"""
struct SharedPhylogeneticInfo{C <: Convention, N} <: TreeSimilarity
    convention::C
    normalize::N
end

SharedPhylogeneticInfo(; convention = TreeDistConvention(), normalize = false) =
    SharedPhylogeneticInfo(Convention(convention), normalize)

function _compare(::SharedPhylogeneticInfo, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    s1, s2 = splits(t1, index), splits(t2, index)
    table = SplitInfoTable(length(index))
    return _matchedphylogeneticinfo(_sharedphylogeneticinfo, s1, s2, table)
end

normalizerinfo(::SharedPhylogeneticInfo, ::Convention, tree) =
    _treephylogeneticinfo(tree)

function normalizer(comparison::SharedPhylogeneticInfo, convention::Convention, t1, t2)
    return (
        normalizerinfo(comparison, convention, t1) +
            normalizerinfo(comparison, convention, t2)
    ) / 2
end

"""
    MatchingSplitInfoDistance(; convention = :treedist, normalize = false)

The matching split information distance of Smith (2020), in bits. For each pair of input
splits, the scorer constructs the two splits supported by taxa whose membership agrees or
differs between the inputs. The more informative of those candidates supplies the pair's
similarity. A linear assignment maximizes the total similarity between non-trivial splits.
The distance is

```math
MSID(T_1, T_2) = I(T_1) + I(T_2) - 2 MSI(T_1, T_2),
```

where `I(T)` is the sum of [`splitinfo`](@ref) over a tree's non-trivial splits and `MSI`
is the optimal matched similarity. An unmatched split contributes its full information to
the distance.

With `normalize = true`, the distance is divided by `I(T₁) + I(T₂)`, following TreeDist.
A function such as `max` or `min` combines the per-tree totals instead, and a number is used
directly as the divisor. Both conventions currently compute the same definition.

Computing all split-pair scores takes `O(s₁ s₂ ceil(n / 64))` time, followed by a dense
linear assignment solve, where `n` is the taxon count and `s₁`, `s₂` are the split counts.

Smith, M.R. (2020). *Information theoretic Generalized Robinson-Foulds metrics for
comparing phylogenetic trees.* Bioinformatics 36(20): 5007–5013.
"""
struct MatchingSplitInfoDistance{C <: Convention, N} <: TreeMetric
    convention::C
    normalize::N
end

MatchingSplitInfoDistance(; convention = TreeDistConvention(), normalize = false) =
    MatchingSplitInfoDistance(Convention(convention), normalize)

function _compare(::MatchingSplitInfoDistance, ::Convention, t1, t2)
    index = taxonindex(t1, t2)
    s1, s2 = splits(t1, index), splits(t2, index)
    table = SplitInfoTable(length(index))
    similarity = _matchedphylogeneticinfo(_matchingsplitinfo, s1, s2, table)
    return _splitwiseinfo(table, s1) + _splitwiseinfo(table, s2) - 2 * similarity
end

normalizerinfo(::MatchingSplitInfoDistance, ::Convention, tree) =
    _treephylogeneticinfo(tree)
