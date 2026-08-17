"""
    TreeMetric <: Distances.SemiMetric

Supertype for distances between two phylogenetic trees.

Subtyping `Distances.SemiMetric` asserts that a metric is symmetric and vanishes on
identical trees, which is what lets `Distances.pairwise` evaluate one triangle and take
the diagonal to be zero. Several tree distances also satisfy the triangle inequality; that
is documented per metric rather than promised here, because not all of them do.

A metric is a value rather than a function: everything that selects a variant of the
computation — the [`Convention`](@ref), whether the result is normalized, and any
metric-specific parameters — is a field, so a constructed metric specifies the computation
completely and needs no keyword arguments when applied.

    RobinsonFoulds(convention = PrimaryConvention(), normalize = true)(t1, t2)

Use [`TreeSimilarity`](@ref) instead for quantities that are largest on identical trees.

# Implementing a metric

A concrete subtype must define

    PhyloDistances._compare(metric, ::TreeDistConvention, t1, t2)

and either carry `convention` and `normalize` fields or override [`convention`](@ref) and
[`isnormalized`](@ref). It may additionally define

  - `_compare(metric, ::PrimaryConvention, t1, t2)` where the primary source differs from
    TreeDist;
  - [`normalizer`](@ref), to support `normalize = true`;
  - [`requiresrooted`](@ref), if the metric is defined on rooted trees;
  - `Distances.result_type`, if the result is not `Float64`.
"""
abstract type TreeMetric <: Distances.SemiMetric end

"""
    TreeSimilarity

Supertype for similarities between two phylogenetic trees: quantities that are *largest*
on identical trees.

These deliberately sit outside the Distances.jl hierarchy. `Distances.PreMetric` requires
`d(x, x) == 0`, which a similarity does not satisfy, and subtyping `SemiMetric` would be
worse than merely inaccurate: `Distances.pairwise` takes the zero diagonal on faith rather
than computing it, so a similarity declared a `SemiMetric` reports zero self-similarity.

`Distances.pairwise` and `Distances.evaluate` are nonetheless defined for similarities, so
they are used exactly like a [`TreeMetric`](@ref); only the type hierarchy differs.

Implementing one follows [`TreeMetric`](@ref).
"""
abstract type TreeSimilarity end

"""
    TreeComparison

Any tree comparison, whether a [`TreeMetric`](@ref) or a [`TreeSimilarity`](@ref). Traits
common to both are defined on this union.
"""
const TreeComparison = Union{TreeMetric,TreeSimilarity}

"""
    Convention

Supertype for the published formulations a metric can be computed under.

Several tree metrics have more than one formulation in the literature, differing in
normalization, in the treatment of trivial splits, or in the base of the logarithm. The
convention selects which of them is computed, so that a value can be traced back to a
specific definition.

Instances are [`TreeDistConvention`](@ref) and [`PrimaryConvention`](@ref). Constructors
taking a convention also accept the shorthand symbols `:treedist` and `:primary`.
"""
abstract type Convention end

"""
    TreeDistConvention()

The formulation used by the R package
[TreeDist](https://github.com/ms609/TreeDist), and the default throughout this package.

Choosing one reference implementation for every metric means values are directly
comparable against it.
"""
struct TreeDistConvention <: Convention end

"""
    PrimaryConvention()

The formulation given by the source that first defined the metric, where that differs
from [`TreeDistConvention`](@ref).

A metric whose two conventions coincide accepts this and returns the same value; a metric
for which only the TreeDist formulation is implemented throws.
"""
struct PrimaryConvention <: Convention end

Base.Symbol(::TreeDistConvention) = :treedist
Base.Symbol(::PrimaryConvention) = :primary

"""
    Convention(name::Symbol) -> Convention

Resolve the shorthand `:treedist` or `:primary` to a [`Convention`](@ref). A `Convention`
is returned unchanged, so a metric constructor may accept either spelling.
"""
function Convention(name::Symbol)
    name === :treedist && return TreeDistConvention()
    name === :primary && return PrimaryConvention()
    throw(ArgumentError(
        "unknown convention $(repr(name)); expected :treedist or :primary"
    ))
end
Convention(convention::Convention) = convention

"""
    convention(comparison) -> Convention

The [`Convention`](@ref) `comparison` computes under.

Reads the `convention` field; override for a metric that stores it differently.
"""
convention(comparison::TreeComparison) = comparison.convention

"""
    isnormalized(comparison) -> Bool

Whether `comparison` divides its result by [`normalizer`](@ref).

Reads the `normalize` field; override for a metric that stores it differently.
"""
isnormalized(comparison::TreeComparison) = comparison.normalize

"""
    requiresrooted(comparison) -> Bool

Whether `comparison` is defined on rooted trees. Defaults to `false`; most tree metrics
are defined on unrooted trees.

The rooting of the inputs is reconciled against this trait before a metric is applied, so
a metric's own implementation may assume it receives trees of the rooting it declares.
"""
requiresrooted(::TreeComparison) = false

"""
    issimilarity(comparison) -> Bool

Whether larger values of `comparison` mean *more* similar trees, i.e. whether it is a
[`TreeSimilarity`](@ref) rather than a [`TreeMetric`](@ref).
"""
issimilarity(::TreeMetric) = false
issimilarity(::TreeSimilarity) = true

"""
    normalizer(comparison, convention, t1, t2) -> Real

The quantity a normalizing comparison divides its result by, typically the largest value
it can take on trees with this taxon set.

Throws unless the metric defines a normalization.
"""
function normalizer(comparison::TreeComparison, ::Convention, t1, t2)
    throw(ArgumentError(
        "no normalization is defined for $(nameof(typeof(comparison))); " *
        "construct it with normalize = false"
    ))
end

"""
    _compare(comparison, convention, t1, t2) -> Real

The unnormalized comparison of `t1` and `t2`. This is the single method a concrete metric
must implement; applying the metric wraps it with normalization.
"""
function _compare(comparison::TreeComparison, convention::Convention, t1, t2)
    throw(ArgumentError(
        "$(nameof(typeof(comparison))) does not implement the " *
        "$(repr(Symbol(convention))) convention"
    ))
end

function _apply(comparison::TreeComparison, t1, t2)
    conv = convention(comparison)
    raw = _compare(comparison, conv, t1, t2)
    isnormalized(comparison) || return raw
    return raw / normalizer(comparison, conv, t1, t2)
end

(metric::TreeMetric)(t1, t2) = _apply(metric, t1, t2)
(similarity::TreeSimilarity)(t1, t2) = _apply(similarity, t1, t2)

Distances.evaluate(similarity::TreeSimilarity, t1, t2) = similarity(t1, t2)

# Distances.jl infers the result type from `one(eltype(a))`, which assumes observations are
# numeric containers. Trees are not, so every comparison must answer for itself.
Distances.result_type(::TreeComparison, ::Type, ::Type) = Float64

"""
    pairwise(comparison, trees)

Compare every pair of `trees`, returning a square matrix `D` with
`D[i, j] == comparison(trees[i], trees[j])`.

For a [`TreeMetric`](@ref) this is `Distances.pairwise`, which evaluates the strict upper
triangle and fills the diagonal with zeros. For a [`TreeSimilarity`](@ref) the diagonal is
computed, since self-similarity is not zero.

The result is always 1-based, even for an offset input; indices refer to positions in
`trees`.
"""
function Distances.pairwise(similarity::TreeSimilarity, trees::AbstractVector)
    return Distances.pairwise(similarity, trees, trees)
end

function Distances.pairwise(
    similarity::TreeSimilarity, a::AbstractVector, b::AbstractVector
)
    T = Distances.result_type(similarity, eltype(a), eltype(b))
    R = Matrix{T}(undef, length(a), length(b))
    ia, ib = collect(eachindex(a)), collect(eachindex(b))
    for j in eachindex(ib), i in eachindex(ia)
        R[i, j] = similarity(a[ia[i]], b[ib[j]])
    end
    return R
end
