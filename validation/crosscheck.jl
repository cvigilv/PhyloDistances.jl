#!/usr/bin/env julia
#
# Checks that this package's values are identical to the R implementations it reproduces,
# not merely close.
#
#     julia --project=validation validation/crosscheck.jl [ncases]
#
# Generates deliberately awkward trees, computes each metric here and in R, and compares
# integers exactly and floating-point values bitwise. Writes validation/report.md.
#
# R is called through RCall, so values cross as machine doubles and integers rather than as
# text. That matters: R's `as.numeric` does not reliably round-trip its own `%.17g` output,
# so a comparison routed through a file has to be written in one particular direction to
# stay honest. Passing the values in memory removes the question.

using Logging: Logging
using PhyloDistances
using PhyloDistances: NewickTree
using Printf
using Random

const HERE = @__DIR__

if Sys.which("R") === nothing
    error("R not found on PATH; the crosscheck compares against TreeDist and Quartet")
end

# RCall records where R lives when it is built, and a Nix rebuild moves that path. Asking
# the R on PATH where it lives keeps the two in step without rebuilding.
haskey(ENV, "R_HOME") || (ENV["R_HOME"] = strip(read(`R RHOME`, String)))
using RCall

star(n) = readnw("(" * join(("T$i:1.0" for i in 1:n), ",") * ");")

"""A caterpillar: the most unbalanced shape on `n` taxa."""
function caterpillar(n)
    nw = "T$n:1.0"
    for i in (n - 1):-1:3
        nw = "($nw,T$i:1.0):1.0"
    end
    return readnw("(T1:1.0,T2:1.0,$nw);")
end

"""Collapse `k` internal branches, turning a binary tree polytomous."""
function collapse!(rng, tree, k)
    for _ in 1:k
        internal = [n for n in NewickTree.prewalk(tree)
                    if !NewickTree.isleaf(n) && n !== tree && isdefined(n, :parent)]
        isempty(internal) && break
        v = rand(rng, internal)
        p = v.parent
        i = findfirst(c -> c === v, p.children)
        i === nothing && continue
        deleteat!(p.children, i)
        for c in v.children
            push!(p, c)
        end
    end
    return tree
end

"""Root an unrooted tree by collapsing its degree-3 root to degree 2."""
function asrooted(tree)
    kids = collect(NewickTree.children(tree))
    length(kids) < 3 && return tree
    root = NewickTree.Node(UInt16(60000), NewickTree.NewickData())
    rest = NewickTree.Node(UInt16(60001), NewickTree.NewickData(1.0, 0.0, ""))
    push!(root, kids[1])
    for kid in kids[2:end]
        push!(rest, kid)
    end
    push!(root, rest)
    return root
end

"""
Tree pairs chosen to stress the cases where implementations usually diverge: no splits at
all, maximal imbalance, polytomies, rooted against unrooted, and identical against
maximally different.
"""
function cases(rng, extra::Int)
    out = Tuple{Any,Any,String}[]
    for n in (4, 5, 6, 7, 8, 12, 20, 50, 120, 400)
        push!(out, (star(n), star(n), "star vs star, n=$n"))
        push!(out, (star(n), caterpillar(n), "star vs caterpillar, n=$n"))
        push!(out, (caterpillar(n), caterpillar(n), "caterpillar identical, n=$n"))
        push!(out, (caterpillar(n), perturb(rng, caterpillar(n), 1), "caterpillar +1 NNI, n=$n"))

        base = randomtree(rng, n)
        push!(out, (base, base, "identical, n=$n"))
        push!(out, (base, perturb(rng, base, 10n), "maximally perturbed, n=$n"))
        push!(out, (asrooted(randomtree(rng, n)), asrooted(randomtree(rng, n)),
                    "rooted vs rooted, n=$n"))
        push!(out, (asrooted(randomtree(rng, n)), randomtree(rng, n),
                    "rooted vs unrooted, n=$n"))

        for frac in (0.2, 0.5, 0.9)
            k = max(1, round(Int, frac * (n - 3)))
            p = collapse!(rng, randomtree(rng, n), k)
            q = collapse!(rng, randomtree(rng, n), k)
            push!(out, (p, q, "polytomies both sides ($frac), n=$n"))
            push!(out, (p, randomtree(rng, n), "polytomies one side ($frac), n=$n"))
        end
    end

    for _ in 1:extra
        n = rand(rng, 4:60)
        a = randomtree(rng, n)
        b = perturb(rng, a, rand(rng, 0:(4n)))
        rand(rng, Bool) && (a = asrooted(a))
        rand(rng, Bool) && (b = asrooted(b))
        rand(rng, Bool) && collapse!(rng, a, rand(rng, 1:max(1, n ÷ 3)))
        push!(out, (a, b, "random, n=$n"))
    end
    return out
end

"""
Load the R side and define the one function the comparison calls.

`QuartetStatus` counts the four-taxon subsets both trees resolve the same way (`s`), resolve
in conflicting ways (`d`), resolve in only one tree (`r1`, `r2`), and leave unresolved in
both (`u`), out of `Q` in total. It also reports `N = 2 * Q`, which is never read here: it
overflows R's 32-bit integers before `Q` does.
"""
function setupR()
    # TreeDist and Quartet both export `RobinsonFoulds`, meaning different things, and
    # whichever is attached second wins. Every call is qualified so load order cannot
    # decide which function runs.
    R"""
    suppressMessages(library(ape))
    suppressMessages(loadNamespace("TreeDist"))
    suppressMessages(loadNamespace("Quartet"))

    compare <- function(nw1, nw2) {
      t1 <- read.tree(text = nw1)
      t2 <- read.tree(text = nw2)
      st <- Quartet::QuartetStatus(c(t1, t2))[2, ]
      c(rf     = TreeDist::RobinsonFoulds(t1, t2),
        rfnorm = TreeDist::RobinsonFoulds(t1, t2, normalize = TRUE),
        Q      = st[["Q"]],
        d      = st[["d"]],
        r1     = st[["r1"]],
        r2     = st[["r2"]])
    }
    """
    return rcopy(R"""
      paste("TreeDist", packageVersion("TreeDist"),
            "| Quartet", packageVersion("Quartet"),
            "| ape", packageVersion("ape"),
            "| R", paste0(R.version$major, ".", R.version$minor))
    """)
end

# Half the cases pass a rooted tree to an unrooted metric on purpose, and each one warns
# that the root position is ignored. That is the behaviour under test, not something to
# report hundreds of times. R's own warnings stay visible: they are raised outside this.
quiet(f) = Logging.with_logger(f, Logging.SimpleLogger(stderr, Logging.Error))

"""
Compare one pair, returning the mismatch descriptions per quantity — empty where the two
implementations agree — and whether the reference normalized Robinson-Foulds was `NaN`.

`NaN` counts as agreement only against `NaN`, which is what a normalized distance gives when
neither tree carries a split and the divisor is zero.
"""
function checkpair(t1, t2, label, nw1, nw2)
    rf, rfnorm, q, d, r1, r2 = rcopy(R"compare($nw1, $nw2)")
    jrf, jrfn, jq, jqn = quiet() do
        (
            RobinsonFoulds()(t1, t2),
            RobinsonFoulds(; normalize = true)(t1, t2),
            QuartetDistance()(t1, t2),
            QuartetDistance(; normalize = true)(t1, t2),
        )
    end
    bad = Pair{Symbol,String}[]

    jrf == rf || push!(bad, :rf => "$label: R=$rf here=$jrf")

    if isnan(rfnorm)
        isnan(jrfn) || push!(bad, :rfnorm => "$label: R=NaN here=$(repr(jrfn))")
    else
        jrfn === rfnorm || push!(bad, :rfnorm => "$label: R=$(repr(rfnorm)) here=$(repr(jrfn))")
    end

    # A quartet that only one tree resolves counts as a difference here, so the reference
    # value is `d + r1 + r2` rather than `d` alone.
    refq = Int(d + r1 + r2)
    jq == refq || push!(bad, :quartet => "$label: R=$refq (d=$d r1=$r1 r2=$r2) here=$jq")

    # Both sides must also agree on how many quartets there are to divide by.
    ntaxa = length(quiet(() -> taxonindex(t1, t2)))
    Int(q) == binomial(ntaxa, 4) ||
        push!(bad, :quartetq => "$label: R Q=$q here=$(binomial(ntaxa, 4))")

    refqn = refq / q
    jqn === refqn || push!(bad, :quartetnorm => "$label: R=$(repr(refqn)) here=$(repr(jqn))")

    return bad, isnan(rfnorm)
end

const QUANTITIES = [
    :rf => ("`RobinsonFoulds()`", "TreeDist", "exact integer"),
    :rfnorm => ("`RobinsonFoulds(normalize = true)`", "TreeDist", "bitwise float"),
    :quartet => ("`QuartetDistance()`", "Quartet", "exact integer"),
    :quartetq => ("quartet count `Q`", "Quartet", "exact integer"),
    :quartetnorm => ("`QuartetDistance(normalize = true)`", "Quartet", "bitwise float"),
]

function main()
    extra = isempty(ARGS) ? 1000 : parse(Int, ARGS[1])
    rng = Xoshiro(20260817)
    pairs = cases(rng, extra)
    @info "generated $(length(pairs)) cases"

    versions = setupR()
    @info "R ready" versions

    mismatches = Dict(first(q) => String[] for q in QUANTITIES)
    nanmatched = 0
    total = 0

    for (t1, t2, label) in pairs
        nw1, nw2 = NewickTree.nwstr(t1), NewickTree.nwstr(t2)
        total += 1
        bad, wasnan = checkpair(t1, t2, label, nw1, nw2)
        wasnan && (nanmatched += 1)
        for (quantity, message) in bad
            push!(mismatches[quantity], message)
        end
        total % 250 == 0 && @info "compared $total / $(length(pairs))"
    end

    io = IOBuffer()
    println(io, "# Agreement with TreeDist and Quartet\n")
    println(io, "Generated by `julia --project=validation validation/crosscheck.jl`.\n")
    println(io, "Integers are compared exactly and floating-point values bitwise — no ")
    println(io, "tolerance anywhere. Values cross from R through RCall as machine numbers, ")
    println(io, "so nothing is routed through text. Cases cover trees with no splits at ")
    println(io, "all, maximal imbalance, polytomies on one and both sides, rooted against ")
    println(io, "unrooted, and identical against maximally different, alongside random ")
    println(io, "trees.\n")
    println(io, "- Julia ", VERSION, ", ", strip(versions))
    println(io, "- Cases compared: ", total)
    println(io)
    println(io, "| quantity | reference | comparison | mismatches |")
    println(io, "|---|---|---|---:|")
    for (quantity, (name, pkg, kind)) in QUANTITIES
        println(io, "| ", name, " | ", pkg, " | ", kind, " | ", length(mismatches[quantity]), " |")
    end
    println(io)
    println(io, "`RobinsonFoulds(normalize = true)` was `NaN` in ", nanmatched,
            " cases, where neither tree carries a split and the divisor is zero.\n")
    println(io, "The quartet distance counts a four-taxon subset that only one tree ")
    println(io, "resolves as a difference, so the reference value is Quartet's `d + r1 + r2`. ")
    println(io, "`Q`, the number of four-taxon subsets, is checked against `binomial(n, 4)` ")
    println(io, "so that both sides agree on the divisor as well as the distance.\n")

    failed = sum(length, values(mismatches))
    if failed == 0
        println(io, "Every value is identical.")
    else
        println(io, "## Mismatches\n")
        for (quantity, _) in QUANTITIES, m in first(mismatches[quantity], 20)
            println(io, "- ", m)
        end
    end
    write(joinpath(HERE, "report.md"), String(take!(io)))

    @info "compared $total cases" mismatches = failed
    failed == 0 || error("values differ from the R reference")
    return nothing
end

main()
