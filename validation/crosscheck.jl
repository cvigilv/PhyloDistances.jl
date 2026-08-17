#!/usr/bin/env julia
#
# Checks that this package's values are identical to TreeDist's, not merely close.
#
#     julia --project=validation validation/crosscheck.jl [ncases]
#
# Generates deliberately awkward trees, computes each metric here and in R, and compares
# integers exactly and floating-point values bitwise. Writes validation/report.md.
#
# Direction matters. R's `as.numeric` does not reliably round-trip its own `%.17g` output —
# it can land half an ulp away — so R writes the values and Julia parses and compares them.
# Comparing the other way measures R's parser rather than either implementation.

using PhyloDistances
using PhyloDistances: NewickTree
using Printf
using Random

const HERE = @__DIR__

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

function main()
    extra = isempty(ARGS) ? 1000 : parse(Int, ARGS[1])
    rng = Xoshiro(20260817)
    pairs = cases(rng, extra)

    jl = joinpath(HERE, "julia_cases.tsv")
    open(jl, "w") do io
        println(io, "label\tnw1\tnw2")
        for (a, b, label) in pairs
            println(io, label, "\t", NewickTree.nwstr(a), "\t", NewickTree.nwstr(b))
        end
    end
    @info "wrote $(length(pairs)) cases"

    if Sys.which("Rscript") === nothing
        @error "Rscript not found; cannot compare against TreeDist"
        return nothing
    end

    rvals = joinpath(HERE, "r_values.tsv")
    info = read(`Rscript $(joinpath(HERE, "treedist_values.R")) $jl $rvals`, String)
    @info "TreeDist finished" info = strip(info)

    rf_bad = String[]
    norm_bad = String[]
    nan_matched = 0
    total = 0

    for row in split(read(rvals, String), '\n'; keepempty = false)[2:end]
        label, nw1, nw2, rrf, rnorm = split(row, '\t')
        total += 1
        t1, t2 = readnw(String(nw1)), readnw(String(nw2))

        jrf = RobinsonFoulds()(t1, t2)
        jrf == parse(Int, rrf) ||
            push!(rf_bad, "$label: TreeDist=$rrf here=$jrf")

        jnm = RobinsonFoulds(; normalize = true)(t1, t2)
        if rnorm == "NaN"
            isnan(jnm) ? (nan_matched += 1) :
                push!(norm_bad, "$label: TreeDist=NaN here=$(repr(jnm))")
        else
            rv = parse(Float64, rnorm)
            jnm === rv ||
                push!(norm_bad, "$label: TreeDist=$(repr(rv)) here=$(repr(jnm))")
        end
    end

    io = IOBuffer()
    println(io, "# Agreement with TreeDist\n")
    println(io, "Generated by `julia --project=validation validation/crosscheck.jl`.\n")
    println(io, "Integers are compared exactly and floating-point values bitwise — no ")
    println(io, "tolerance anywhere. Cases cover trees with no splits at all, maximal ")
    println(io, "imbalance, polytomies on one and both sides, rooted against unrooted, and ")
    println(io, "identical against maximally different, alongside random trees.\n")
    println(io, "- Julia ", VERSION, ", ", strip(info))
    println(io, "- Cases compared: ", total)
    println(io)
    println(io, "| quantity | comparison | mismatches |")
    println(io, "|---|---|---:|")
    println(io, "| `RobinsonFoulds()` | exact integer | ", length(rf_bad), " |")
    println(io, "| `RobinsonFoulds(normalize = true)` | bitwise float | ", length(norm_bad), " |")
    println(io)
    println(io, "`NaN` agreed with `NaN` in ", nan_matched,
            " cases, where neither tree carries a split and the normalizer is zero.\n")
    if isempty(rf_bad) && isempty(norm_bad)
        println(io, "Every value is identical.")
    else
        println(io, "## Mismatches\n")
        for m in first(vcat(rf_bad, norm_bad), 40)
            println(io, "- ", m)
        end
    end
    write(joinpath(HERE, "report.md"), String(take!(io)))

    @info "compared $total cases" rf_mismatches = length(rf_bad) norm_mismatches = length(norm_bad)
    isempty(rf_bad) && isempty(norm_bad) || error("values differ from TreeDist")
    return nothing
end

main()
