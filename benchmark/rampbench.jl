#!/usr/bin/env julia
#
# A doubling taxon ramp for the metrics that build a score matrix over the taxa, run
# against the same R references as run.jl.
#
#     julia --project=benchmark benchmark/rampbench.jl
#
# Writes the trees both sides read, runs each, and renders ramp.md. R is optional: without
# a working Rscript the Julia timings are reported alone.
#
# The sizes are powers of two on purpose, which is what makes this worth running alongside
# run.jl's round numbers rather than folding into it. `QuartetDistance` holds an n x n
# table indexed at scattered [x, y]; a leading dimension that is a power of two sends whole
# families of those accesses to one cache set, and the resulting cliff is invisible at
# 10/50/200/1000 taxa. Anything else laid out as a square table is liable to the same
# thing, so a ramp that lands squarely on the bad sizes is the one that finds it.
#
# The quartet distance at 1024 taxa takes seconds; a full run is not quick.

using BenchmarkTools
using PhyloDistances
using PhyloDistances: NewickTree
using Printf
using Random

const HERE = @__DIR__
const TREEDIR = joinpath(HERE, "ramptrees")
const RFILE = joinpath(HERE, "results_ramp_r.tsv")
const JULIAFILE = joinpath(HERE, "results_ramp_julia.tsv")

const SIZES = (16, 64, 256, 1024)

# Quartet represents quartet counts in 32-bit integers and refuses trees above 477 tips.
const QUARTET_CEILING = 477

# Above this a single quartet call takes seconds, so it is measured once with `@timed`
# rather than sampled: the run-to-run variation is far smaller than the quantity itself.
const LONGCALL = 512

"""Write the Newick files both implementations read, so neither gets different inputs."""
function writetrees()
    mkpath(TREEDIR)
    rng = Xoshiro(20260821)

    for n in SIZES
        a = randomtree(rng, n)
        b = perturb(rng, a, n ÷ 4)
        write(joinpath(TREEDIR, "pair_$(n)_a.nwk"), NewickTree.nwstr(a) * "\n")
        write(joinpath(TREEDIR, "pair_$(n)_b.nwk"), NewickTree.nwstr(b) * "\n")
    end
    return nothing
end

readpair(n) = (
    readnw(read(joinpath(TREEDIR, "pair_$(n)_a.nwk"), String)),
    readnw(read(joinpath(TREEDIR, "pair_$(n)_b.nwk"), String)),
)

"""
Time one comparison, recording the value alongside so the two implementations can be
checked against each other on exactly the trees benchmarked. `allocs` is `nothing` for a
call timed once rather than sampled, which is how `@timed` reports it.
"""
function measure(metric, comparison, n, t1, t2)
    if metric == "quartet" && n >= LONGCALL
        comparison(t1, t2)                      # warm up; discard
        stat = @timed comparison(t1, t2)
        return (; metric, n, seconds = stat.time, megabytes = stat.bytes / 1024^2,
                allocs = nothing, value = Float64(stat.value))
    end

    trial = @benchmark $comparison($t1, $t2) samples = 100 seconds = 10
    return (; metric, n, seconds = minimum(trial).time / 1e9,
            megabytes = minimum(trial).memory / 1024^2,
            allocs = minimum(trial).allocs, value = Float64(comparison(t1, t2)))
end

"""
Time Jaccard-Robinson-Foulds and Nye similarity at their shared `k = 1,
allowconflict = true` scoring, and the quartet distance, on every size of the ramp.
"""
function benchmarkjulia()
    rows = NamedTuple[]

    for n in SIZES
        t1, t2 = readpair(n)
        for (metric, comparison) in (
            "jrf" => JaccardRobinsonFoulds(),
            "nye" => NyeSimilarity(),
            "quartet" => QuartetDistance(),
        )
            @info "$metric at $n taxa"
            push!(rows, measure(metric, comparison, n, t1, t2))
        end
    end

    return rows
end

"""Write the Julia side in the same shape the R script writes its own."""
function writetsv(path, rows)
    open(path, "w") do io
        println(io, "metric\tn\tseconds\tmegabytes\tallocs\tvalue")
        for r in rows
            @printf(io, "%s\t%d\t%.9g\t%.6g\t%s\t%.12g\n", r.metric, r.n, r.seconds,
                    r.megabytes, isnothing(r.allocs) ? "NA" : string(r.allocs), r.value)
        end
    end
    return nothing
end

function readtsv(path)
    isfile(path) || return NamedTuple[]
    lines = readlines(path)
    length(lines) > 1 || return NamedTuple[]
    header = split(strip(lines[1]), '\t')
    return map(lines[2:end]) do line
        vals = Dict(k => v for (k, v) in zip(header, split(strip(line), '\t')))
        (
            metric = strip(vals["metric"], '"'),
            n = parse(Int, vals["n"]),
            seconds = parse(Float64, vals["seconds"]),
            megabytes = parse(Float64, vals["megabytes"]),
            value = parse(Float64, vals["value"]),
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
    return ours <= theirs ? @sprintf("%.2f× faster", theirs / ours) :
           @sprintf("%.2f× slower", ours / theirs)
end

const TITLES = Dict(
    "jrf" => "Jaccard-Robinson-Foulds",
    "nye" => "Nye similarity",
    "quartet" => "Quartet distance",
)

const REFERENCES = Dict(
    "jrf" => "TreeDist", "nye" => "TreeDist", "quartet" => "Quartet",
)

function section(io, metric, rows, rrows)
    println(io, "## ", TITLES[metric], "\n")
    println(io, "| taxa | PhyloDistances | ", REFERENCES[metric],
            " | ratio | Julia alloc | Julia allocs | agree |")
    println(io, "|-----:|---------------:|---------:|------:|------------:|-------------:|:------|")

    for row in filter(x -> x.metric == metric, rows)
        match = findfirst(x -> x.metric == metric && x.n == row.n, rrows)
        rt = isnothing(match) ? nothing : rrows[match].seconds
        agree = if isnothing(match)
            "—"
        elseif isapprox(rrows[match].value, row.value; atol = 1e-9, rtol = 1e-6)
            "yes"
        else
            @sprintf("**NO** (R %.6f vs %.6f)", rrows[match].value, row.value)
        end
        @printf(io, "| %d | %s | %s | %s | %.2f MB | %s | %s |\n",
            row.n, prettytime(row.seconds), isnothing(rt) ? "—" : prettytime(rt),
            prettyratio(row.seconds, rt), row.megabytes,
            isnothing(row.allocs) ? "—" : string(row.allocs), agree)
    end
    println(io)
    return nothing
end

function report(rows, rrows, rinfo)
    io = IOBuffer()
    println(io, "# PhyloDistances.jl on a doubling taxon ramp\n")
    println(io, "Generated by `julia --project=benchmark benchmark/rampbench.jl`. ")
    println(io, "Companion to `results.md`, which uses round taxon counts; this one uses ")
    println(io, "powers of two, where a table whose leading dimension is a power of two ")
    println(io, "collides in the cache and the round sizes show nothing.\n")
    println(io, "Both implementations read the same Newick files from `benchmark/ramptrees/`, ")
    println(io, "and neither timing includes parsing — only the distance computation.\n")
    println(io, "**Memory is measured differently on each side**, so only the Julia column ")
    println(io, "is reported here; see `results.md` on why the two are not comparable.\n")
    println(io, "Quartet refuses trees above ", QUARTET_CEILING, " tips, where the number of ")
    println(io, "four-taxon subsets outgrows the 32-bit integers it counts them in, so the ")
    println(io, "1024-taxon row has no reference to compare against.\n")

    println(io, "## Environment\n")
    println(io, "- Julia ", VERSION, " on ", Sys.MACHINE)
    for line in unique(strip.(rinfo))
        println(io, "- ", line)
    end
    println(io)

    for metric in ("jrf", "nye", "quartet")
        section(io, metric, rows, rrows)
    end

    println(io, "`results_ramp_julia.tsv` and `results_ramp_r.tsv` hold the same numbers ")
    println(io, "row by row, and are regenerated on every run. Keep a copy of one before a ")
    println(io, "change to measure the change against; the trees are written from a fixed ")
    println(io, "seed, so two runs compare like for like.")

    write(joinpath(HERE, "ramp.md"), String(take!(io)))
    return nothing
end

"""Run the R benchmark script, returning its version banner or `nothing` if it failed."""
function runR()
    try
        return read(
            `Rscript $(joinpath(HERE, "rampbench.R")) $TREEDIR $RFILE $(join(SIZES, ",")) $QUARTET_CEILING`,
            String,
        )
    catch err
        @warn "rampbench.R failed; reporting Julia timings only" err
        return nothing
    end
end

function main()
    @info "Writing trees to $TREEDIR"
    writetrees()

    @info "Benchmarking PhyloDistances"
    rows = benchmarkjulia()

    banner = if Sys.which("Rscript") === nothing
        @warn "Rscript not found; reporting Julia timings only"
        nothing
    else
        @info "Benchmarking TreeDist and Quartet"
        runR()
    end

    writetsv(JULIAFILE, rows)
    report(rows, readtsv(RFILE), isnothing(banner) ? String[] : [banner])
    @info "Wrote $(joinpath(HERE, "ramp.md"))"
    return nothing
end

main()
