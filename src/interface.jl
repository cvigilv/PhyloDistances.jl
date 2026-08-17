"""
    TreeMetric

Supertype for comparisons between two phylogenetic trees.

A metric is a value rather than a function: parameters that select a variant of the
computation — exponents, interpolation weights, tie-breaking rules — are fields of the
concrete type, so a constructed metric specifies the computation completely. Apply it with
[`compare`](@ref).

Every metric is symmetric: `compare(metric, t1, t2) == compare(metric, t2, t1)`.
[`pairwise`](@ref) relies on this to fill only one triangle.

# Implementing a metric

A concrete subtype must define

    PhyloDistances._compare(metric, ::TreeDistConvention, t1, t2)

and may additionally define

  - `_compare(metric, ::PrimaryConvention, t1, t2)` where the primary source differs from
    TreeDist;
  - [`normalizer`](@ref), to support `normalize = true`;
  - [`requiresrooted`](@ref), if the metric is defined on rooted trees;
  - [`issimilarity`](@ref), if larger values mean more similar trees.
"""
abstract type TreeMetric end

"""
    Convention

Supertype for the published formulations a metric can be computed under.

Several tree metrics have more than one formulation in the literature, differing in
normalization, in the treatment of trivial splits, or in the base of the logarithm. The
convention selects which of them is computed, so that a value can be traced back to a
specific definition.

Instances are [`TreeDistConvention`](@ref) and [`PrimaryConvention`](@ref). Functions
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
is returned unchanged, so functions may accept either spelling.
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
    requiresrooted(metric) -> Bool

Whether `metric` is defined on rooted trees. Defaults to `false`; most tree metrics are
defined on unrooted trees.

The rooting of the inputs is reconciled against this trait before a metric is applied, so
a metric's own implementation may assume it receives trees of the rooting it declares.
"""
requiresrooted(::TreeMetric) = false

"""
    issimilarity(metric) -> Bool

Whether larger values of `metric` mean *more* similar trees. Defaults to `false`, i.e.
the value is a distance and identical trees score lowest.

This is what distinguishes, say, shared phylogenetic information from the distance
derived from it; it carries no claim that the value satisfies the metric axioms.
"""
issimilarity(::TreeMetric) = false

"""
    normalizer(metric, convention, t1, t2) -> Real

The quantity that [`compare`](@ref) divides by when called with `normalize = true`,
typically the largest value `metric` can take on trees with this taxon set.

Throws unless the metric defines a normalization.
"""
function normalizer(metric::TreeMetric, ::Convention, t1, t2)
    throw(ArgumentError(
        "no normalization is defined for $(nameof(typeof(metric))); " *
        "call with normalize = false"
    ))
end

"""
    _compare(metric, convention, t1, t2) -> Real

The unnormalized comparison of `t1` and `t2`. This is the single method a concrete metric
must implement; [`compare`](@ref) wraps it with convention resolution and normalization.
"""
function _compare(metric::TreeMetric, convention::Convention, t1, t2)
    throw(ArgumentError(
        "$(nameof(typeof(metric))) does not implement the " *
        "$(repr(Symbol(convention))) convention"
    ))
end

"""
    compare(metric, t1, t2; convention = :treedist, normalize = false)

Compare the phylogenetic trees `t1` and `t2` under `metric`.

`convention` selects which published formulation to compute, as either a
[`Convention`](@ref) or the shorthand `:treedist` / `:primary`. With `normalize = true`
the result is divided by [`normalizer`](@ref), placing it on a scale comparable across
taxon set sizes.

Whether a larger value means more or less similar trees depends on the metric; see
[`issimilarity`](@ref).

See [`pairwise`](@ref) to compare every pair in a collection.
"""
function compare(
    metric::TreeMetric, t1, t2; convention = :treedist, normalize::Bool = false
)
    conv = Convention(convention)
    raw = _compare(metric, conv, t1, t2)
    normalize || return raw
    return raw / normalizer(metric, conv, t1, t2)
end

"""
    pairwise(metric, trees; convention = :treedist, normalize = false)

Compare every pair of `trees` under `metric`, returning a square matrix `D` with
`D[i, j] == compare(metric, trees[i], trees[j])`.

The result shares the axes of `trees` along both dimensions. Because metrics are
symmetric only one triangle is evaluated and the other is mirrored, so `metric` is applied
`n(n+1)/2` times for `n` trees. An empty collection yields a 0×0 `Float64` matrix.

Keyword arguments are forwarded to [`compare`](@ref).
"""
function pairwise(metric::TreeMetric, trees::AbstractVector; kwargs...)
    ax = axes(trees, 1)
    isempty(trees) && return similar(trees, Float64, (ax, ax))

    # The diagonal doubles as element-type inference: it is part of the result, so no
    # comparison is computed only to be discarded.
    diagonal = map(t -> compare(metric, t, t; kwargs...), trees)
    D = similar(trees, eltype(diagonal), (ax, ax))
    for i in ax
        D[i, i] = diagonal[i]
    end
    for i in ax, j in ax
        j <= i && continue
        D[i, j] = D[j, i] = compare(metric, trees[i], trees[j]; kwargs...)
    end
    return D
end
