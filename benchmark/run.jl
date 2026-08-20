#!/usr/bin/env julia
#
# Compares this package against the R implementations it reproduces, on identical inputs.
#
#     julia --project=benchmark benchmark/run.jl
#
# Writes the trees both sides read, runs each, and renders results.md. R is optional:
# without a working Rscript the Julia timings are reported alone.
#
# R is called as a subprocess rather than through RCall so that each timing loop runs inside
# R's own clock. The round trip through RCall costs 12-30 µs, which is comparable to the
# fastest quantities measured here and would dominate them.
#
# The quartet distance at 1000 taxa takes minutes; a full run is not quick.

using BenchmarkTools
using PhyloDistances
using PhyloDistances: NewickTree
using Printf
using Random

const HERE = @__DIR__
const TREEDIR = joinpath(HERE, "trees")

const RF_SIZES = (10, 50, 200, 1000)

# Quartet represents quartet counts in 32-bit integers and refuses trees above 477 tips, so
# that is where the comparison stops and the Julia-only sizes take over.
const QUARTET_SIZES = (10, 50, 200, 477)
const QUARTET_LARGE = (700, 1000, 1500)
const QUARTET_CEILING = 477

const PAIR_SIZES = Tuple(sort(unique((RF_SIZES..., QUARTET_SIZES..., QUARTET_LARGE...))))
const COLLECTION_SIZE = 60      # taxa per tree in the all-pairs case
const COLLECTION_COUNT = 40     # trees compared against each other

"""Write the Newick files both implementations read, so neither gets different inputs."""
function writetrees()
    mkpath(TREEDIR)
    rng = Xoshiro(20260817)

    for n in PAIR_SIZES
        a = randomtree(rng, n)
        b = perturb(rng, a, n ÷ 4)
        write(joinpath(TREEDIR, "pair_$(n)_a.nwk"), NewickTree.nwstr(a) * "\n")
        write(joinpath(TREEDIR, "pair_$(n)_b.nwk"), NewickTree.nwstr(b) * "\n")
    end

    base = randomtree(rng, COLLECTION_SIZE)
    open(joinpath(TREEDIR, "collection.nwk"), "w") do io
        for _ in 1:COLLECTION_COUNT
            println(io, NewickTree.nwstr(perturb(rng, base, rand(rng, 1:20))))
        end
    end
    return nothing
end

readpair(n) = (
    readnw(read(joinpath(TREEDIR, "pair_$(n)_a.nwk"), String)),
    readnw(read(joinpath(TREEDIR, "pair_$(n)_b.nwk"), String)),
)

"""Time Robinson-Foulds, reading the same files rather than reusing in-memory trees."""
function benchmarkrf()
    rows = NamedTuple[]

    for n in RF_SIZES
        t1, t2 = readpair(n)
        trial = @benchmark RobinsonFoulds()($t1, $t2) samples = 200 seconds = 3
        push!(rows, (
            case = "pair", n = n, ntrees = 2,
            seconds = minimum(trial).time / 1e9,
            megabytes = minimum(trial).memory / 1024^2,
            allocs = minimum(trial).allocs,
        ))
    end

    trees = [readnw(line) for line in eachline(joinpath(TREEDIR, "collection.nwk")) if !isempty(line)]
    trial = @benchmark pairwise(RobinsonFoulds(), $trees) samples = 50 seconds = 10
    push!(rows, (
        case = "allpairs", n = COLLECTION_SIZE, ntrees = length(trees),
        seconds = minimum(trial).time / 1e9,
        megabytes = minimum(trial).memory / 1024^2,
        allocs = minimum(trial).allocs,
    ))

    return rows
end

"""
Time Jaccard-Robinson-Foulds at its `k = 1, allowconflict = true` default — the same
scoring [`NyeSimilarity`](@ref) uses — recording the value alongside so the two
implementations can be checked against each other on exactly the trees benchmarked.
"""
function benchmarkjrf()
    rows = NamedTuple[]

    for n in RF_SIZES
        t1, t2 = readpair(n)
        trial = @benchmark JaccardRobinsonFoulds()($t1, $t2) samples = 200 seconds = 3
        push!(rows, (
            case = "pair", n = n, ntrees = 2,
            seconds = minimum(trial).time / 1e9,
            megabytes = minimum(trial).memory / 1024^2,
            allocs = minimum(trial).allocs,
            value = JaccardRobinsonFoulds()(t1, t2),
        ))
    end

    return rows
end

"""
Time InfoRobinsonFoulds, recording the value alongside so the two implementations can be
checked against each other on exactly the trees benchmarked.
"""
function benchmarkinforf()
    rows = NamedTuple[]

    for n in RF_SIZES
        t1, t2 = readpair(n)
        trial = @benchmark InfoRobinsonFoulds()($t1, $t2) samples = 200 seconds = 3
        push!(rows, (
            case = "pair", n = n, ntrees = 2,
            seconds = minimum(trial).time / 1e9,
            megabytes = minimum(trial).memory / 1024^2,
            allocs = minimum(trial).allocs,
            value = InfoRobinsonFoulds()(t1, t2),
        ))
    end

    return rows
end

"""
Time the quartet distance, recording the value alongside so the two implementations can be
checked against each other on exactly the trees that were benchmarked.

Above a few hundred taxa a single call takes seconds, so it is measured once with `@timed`
rather than sampled: the run-to-run variation is far smaller than the quantity itself, and
repeating it would multiply an already long benchmark.
"""
function benchmarkquartet()
    rows = NamedTuple[]

    for n in (QUARTET_SIZES..., QUARTET_LARGE...)
        t1, t2 = readpair(n)
        @info "quartet distance at $n taxa ($(binomial(n, 4)) quartets)"

        if n >= 400
            stat = @timed QuartetDistance()(t1, t2)
            push!(rows, (
                case = "quartet", n = n, ntrees = 2,
                seconds = stat.time, megabytes = stat.bytes / 1024^2, value = stat.value,
            ))
        else
            trial = @benchmark QuartetDistance()($t1, $t2) samples = 50 seconds = 5
            push!(rows, (
                case = "quartet", n = n, ntrees = 2,
                seconds = minimum(trial).time / 1e9,
                megabytes = minimum(trial).memory / 1024^2,
                value = QuartetDistance()(t1, t2),
            ))
        end
    end

    return rows
end

function readtsv(path)
    isfile(path) || return NamedTuple[]
    lines = readlines(path)
    length(lines) > 1 || return NamedTuple[]
    header = split(strip(lines[1]), '\t')
    return map(lines[2:end]) do line
        fields = split(strip(line), '\t')
        vals = Dict(k => v for (k, v) in zip(header, fields))
        (
            case = strip(vals["case"], '"'),
            n = parse(Int, vals["n"]),
            ntrees = parse(Int, vals["ntrees"]),
            seconds = parse(Float64, vals["seconds"]),
            megabytes = parse(Float64, vals["megabytes"]),
            value = haskey(vals, "value") ? parse(Float64, vals["value"]) : NaN,
        )
    end
end

prettytime(s) = s < 1e-6 ? @sprintf("%.0f ns", s * 1e9) :
                s < 1e-3 ? @sprintf("%.1f µs", s * 1e6) :
                s < 1.0 ? @sprintf("%.2f ms", s * 1e3) : @sprintf("%.2f s", s)

"""
Express one timing against the other in the direction that reads correctly, so that a
faster result never appears as a fraction of a slower one.
"""
function prettyratio(ours, theirs)
    isnothing(theirs) && return "—"
    return ours <= theirs ? @sprintf("%.1f× faster", theirs / ours) :
           @sprintf("%.1f× slower", ours / theirs)
end

"""Group digits so that quartet counts in the billions stay readable."""
function prettycount(n::Integer)
    s = string(n)
    parts = String[]
    while length(s) > 3
        pushfirst!(parts, s[(end - 2):end])
        s = s[1:(end - 3)]
    end
    pushfirst!(parts, s)
    return join(parts, ",")
end

function report(rf, quartet, jrf, inforf, rrf, rquartet, rjrf, rinforf, rinfo)
    io = IOBuffer()
    println(io, "# PhyloDistances.jl against its R references\n")
    println(io, "Generated by `julia --project=benchmark benchmark/run.jl`.\n")
    println(io, "Both implementations read the same Newick files from `benchmark/trees/`, ")
    println(io, "and neither timing includes parsing — only the distance computation.\n")
    println(io, "- Julia: minimum of a `BenchmarkTools` trial, or a single `@timed` call ")
    println(io, "  where one call already takes seconds.")
    println(io, "- R: median of repeated calls, looped until the clock resolution stops mattering.\n")
    println(io, "**Memory is measured differently on each side and the two columns are not ")
    println(io, "comparable.** Julia reports bytes allocated by the call; R reports the peak ")
    println(io, "its garbage collector saw, which includes everything already resident.\n")

    println(io, "## Environment\n")
    println(io, "- Julia ", VERSION, " on ", Sys.MACHINE)
    for line in unique(strip.(rinfo))
        println(io, "- ", line)
    end
    println(io)

    println(io, "## Robinson-Foulds, one pair of trees\n")
    println(io, "| taxa | PhyloDistances | TreeDist | ratio | Julia alloc | Julia allocs |")
    println(io, "|-----:|---------------:|---------:|------:|------------:|-------------:|")
    for row in filter(x -> x.case == "pair", rf)
        match = findfirst(x -> x.case == "pair" && x.n == row.n, rrf)
        rt = isnothing(match) ? nothing : rrf[match].seconds
        @printf(io, "| %d | %s | %s | %s | %.2f MB | %d |\n",
            row.n, prettytime(row.seconds),
            isnothing(rt) ? "—" : prettytime(rt), prettyratio(row.seconds, rt),
            row.megabytes, row.allocs)
    end

    println(io, "\n## Robinson-Foulds, all pairs within a collection\n")
    for row in filter(x -> x.case == "allpairs", rf)
        match = findfirst(x -> x.case == "allpairs", rrf)
        rt = isnothing(match) ? nothing : rrf[match].seconds
        npairs = row.ntrees * (row.ntrees - 1) ÷ 2
        println(io, "$(row.ntrees) trees of $(row.n) taxa — $npairs pairs\n")
        println(io, "| | time | per pair |")
        println(io, "|---|-----:|---------:|")
        @printf(io, "| PhyloDistances | %s | %s |\n",
            prettytime(row.seconds), prettytime(row.seconds / npairs))
        if !isnothing(rt)
            @printf(io, "| TreeDist | %s | %s |\n", prettytime(rt), prettytime(rt / npairs))
            @printf(io, "| ratio | %s | |\n", prettyratio(row.seconds, rt))
        end
    end

    println(io, "\n## Quartet distance, one pair of trees\n")
    println(io, "Quartet refuses trees above ", QUARTET_CEILING, " tips, where the number of ")
    println(io, "four-taxon subsets outgrows the 32-bit integers it counts them in, so the ")
    println(io, "rows past that size have no reference to compare against.\n")
    println(io, "| taxa | quartets | PhyloDistances | Quartet | ratio | Julia alloc | agree |")
    println(io, "|-----:|---------:|---------------:|--------:|------:|------------:|:------|")
    for row in quartet
        match = findfirst(x -> x.case == "quartet" && x.n == row.n, rquartet)
        rt = isnothing(match) ? nothing : rquartet[match].seconds
        agree = if isnothing(match)
            "—"
        elseif rquartet[match].value == row.value
            "yes"
        else
            @sprintf("**NO** (R %.0f vs %d)", rquartet[match].value, row.value)
        end
        @printf(io, "| %d | %s | %s | %s | %s | %.2f MB | %s |\n",
            row.n, prettycount(binomial(row.n, 4)), prettytime(row.seconds),
            isnothing(rt) ? "—" : prettytime(rt), prettyratio(row.seconds, rt),
            row.megabytes, agree)
    end

    println(io, "\n`QuartetDistance`'s default `:fast` algorithm counts concordant quartets ")
    println(io, "in `O(n³)` without enumerating them (`algorithm = :naive` keeps the exact ")
    println(io, "`O(n⁴)` enumeration as a correctness oracle). Quartet wraps tqDist, which ")
    println(io, "counts the same quantity in `O(n log n)` without enumerating it either — ")
    println(io, "so the remaining gap past 477 taxa is an algorithm gap, not a language one, ")
    println(io, "and it would still widen without bound past where this benchmark stops.")

    println(io, "\n## Jaccard-Robinson-Foulds, one pair of trees\n")
    println(io, "`JaccardRobinsonFoulds()` at its `k = 1, allowconflict = true` default — the ")
    println(io, "scoring [`NyeSimilarity`](@ref) also uses. Both implementations solve an ")
    println(io, "assignment problem over the trees' splits; `agree` uses a tolerance rather ")
    println(io, "than exact equality because the two solvers round differently internally ")
    println(io, "(see `validation/crosscheck.jl`'s `_closeenough`), not because either is wrong.\n")
    println(io, "| taxa | PhyloDistances | TreeDist | ratio | Julia alloc | Julia allocs | agree |")
    println(io, "|-----:|---------------:|---------:|------:|------------:|-------------:|:------|")
    for row in filter(x -> x.case == "pair", jrf)
        match = findfirst(x -> x.case == "pair" && x.n == row.n, rjrf)
        rt = isnothing(match) ? nothing : rjrf[match].seconds
        agree = if isnothing(match)
            "—"
        elseif isapprox(rjrf[match].value, row.value; atol = 1e-9, rtol = 1e-6)
            "yes"
        else
            @sprintf("**NO** (R %.6f vs %.6f)", rjrf[match].value, row.value)
        end
        @printf(io, "| %d | %s | %s | %s | %.2f MB | %d | %s |\n",
            row.n, prettytime(row.seconds),
            isnothing(rt) ? "—" : prettytime(rt), prettyratio(row.seconds, rt),
            row.megabytes, row.allocs, agree)
    end

    println(io, "\n## Info-Robinson-Foulds, one pair of trees\n")
    println(io, "`InfoRobinsonFoulds()`, which weights each split by its phylogenetic ")
    println(io, "information content rather than counting it as one. Unlike ")
    println(io, "Jaccard-Robinson-Foulds, this does not go through an assignment solve — ")
    println(io, "splits are matched by identity, the same relationship classic ")
    println(io, "Robinson-Foulds computes — so `agree` uses a tolerance only because ")
    println(io, "TreeDist floors its result near zero (`.FloorNumericalNoise`), not because ")
    println(io, "either side rounds an optimization differently.\n")
    println(io, "| taxa | PhyloDistances | TreeDist | ratio | Julia alloc | Julia allocs | agree |")
    println(io, "|-----:|---------------:|---------:|------:|------------:|-------------:|:------|")
    for row in filter(x -> x.case == "pair", inforf)
        match = findfirst(x -> x.case == "pair" && x.n == row.n, rinforf)
        rt = isnothing(match) ? nothing : rinforf[match].seconds
        agree = if isnothing(match)
            "—"
        elseif isapprox(rinforf[match].value, row.value; atol = 1e-9, rtol = 1e-6)
            "yes"
        else
            @sprintf("**NO** (R %.6f vs %.6f)", rinforf[match].value, row.value)
        end
        @printf(io, "| %d | %s | %s | %s | %.2f MB | %d | %s |\n",
            row.n, prettytime(row.seconds),
            isnothing(rt) ? "—" : prettytime(rt), prettyratio(row.seconds, rt),
            row.megabytes, row.allocs, agree)
    end

    write(joinpath(HERE, "results.md"), String(take!(io)))
    return nothing
end

"""Run one R benchmark script, returning its version banner or `nothing` if it failed."""
function runR(script, outfile, args...)
    try
        return read(`Rscript $(joinpath(HERE, script)) $TREEDIR $outfile $(args...)`, String)
    catch err
        @warn "$script failed; reporting Julia timings only" err
        return nothing
    end
end

function main()
    @info "Writing trees to $TREEDIR"
    writetrees()

    @info "Benchmarking Robinson-Foulds"
    rf = benchmarkrf()
    @info "Benchmarking the quartet distance"
    quartet = benchmarkquartet()
    @info "Benchmarking Jaccard-Robinson-Foulds"
    jrf = benchmarkjrf()
    @info "Benchmarking InfoRobinsonFoulds"
    inforf = benchmarkinforf()

    rffile = joinpath(HERE, "results_rf_r.tsv")
    qfile = joinpath(HERE, "results_quartet_r.tsv")
    jrffile = joinpath(HERE, "results_jrf_r.tsv")
    inforffile = joinpath(HERE, "results_inforf_r.tsv")
    banners = String[]

    if Sys.which("Rscript") === nothing
        @warn "Rscript not found; reporting Julia timings only"
    else
        @info "Benchmarking TreeDist"
        rfinfo = runR("treedist.R", rffile)
        isnothing(rfinfo) || push!(banners, rfinfo)
        @info "Benchmarking Quartet"
        qinfo = runR("quartet.R", qfile, join(QUARTET_SIZES, ","))
        isnothing(qinfo) || push!(banners, qinfo)
        @info "Benchmarking TreeDist's JaccardRobinsonFoulds"
        jrfinfo = runR("jrf.R", jrffile)
        isnothing(jrfinfo) || push!(banners, jrfinfo)
        @info "Benchmarking TreeDist's InfoRobinsonFoulds"
        inforfinfo = runR("inforf.R", inforffile)
        isnothing(inforfinfo) || push!(banners, inforfinfo)
    end

    report(
        rf, quartet, jrf, inforf,
        readtsv(rffile), readtsv(qfile), readtsv(jrffile), readtsv(inforffile),
        banners,
    )
    @info "Wrote $(joinpath(HERE, "results.md"))"
    return nothing
end

main()
