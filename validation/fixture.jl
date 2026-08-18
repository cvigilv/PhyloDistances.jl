#!/usr/bin/env julia
#
# Sources the committed test fixture from the R implementations this package reproduces.
#
#     julia --project=validation validation/fixture.jl           # check
#     julia --project=validation validation/fixture.jl --write    # regenerate
#
# `test/fixtures/rf_quartet.tsv` is read by the portable test suite, which must pass with
# no R installed. Its expected values are nonetheless not invented here: this script
# computes each one with TreeDist and Quartet and either verifies the committed file
# against them or rewrites it. Checking exits non-zero on any difference, so the fixture
# cannot drift away from the reference unnoticed.
#
# Only the numeric columns are generated. The `provenance` column is written by hand and
# carried through unchanged: it records the independent derivation of each row, which is
# what makes a disagreement between this script and the file informative rather than
# circular.

using Logging: Logging
using PhyloDistances

const HERE = @__DIR__
include(joinpath(HERE, "..", "test", "fixtures", "read.jl"))

if Sys.which("R") === nothing
    error("R not found on PATH; the fixture is generated from TreeDist and Quartet")
end

# RCall records where R lives when it is built, and moving or replacing that installation
# leaves the two out of step. Asking the R on PATH where it lives keeps them together.
haskey(ENV, "R_HOME") || (ENV["R_HOME"] = strip(read(`R RHOME`, String)))
using RCall

"""
Load the R side and define the one function this script calls.

`QuartetStatus` reports the four-taxon subsets the two trees resolve alike (`s`), resolve
in conflicting ways (`d`), resolve in only one of them (`r1`, `r2`), and leave unresolved
in both (`u`), out of `Q` in total. `N = 2 * Q` is never read: it overflows R's 32-bit
integers before `Q` does.

TreeDist and Quartet both export `RobinsonFoulds`, meaning different things, and whichever
is attached second wins. Every call is qualified so load order cannot decide which
function runs.
"""
function setupR()
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

"""
    referencevalues(newick1, newick2) -> NamedTuple

The four expected values for one pair, as the R references give them.

A quartet resolved by only one tree counts as a difference, so the quartet distance is
`d + r1 + r2` rather than `d` alone, and its divisor is `Q`.
"""
function referencevalues(newick1, newick2)
    rf, rfnorm, Q, d, r1, r2 = rcopy(R"compare($newick1, $newick2)")
    quartet = Int(d + r1 + r2)
    return (
        rf = Int(rf),
        rf_normalized = Float64(rfnorm),
        quartet = quartet,
        quartet_normalized = quartet / Q,
    )
end

# `NaN` agrees only with `NaN`, which is what a normalized distance gives when neither tree
# carries a split and the divisor is zero. Everything else is compared bitwise: these are
# the same IEEE operations on both sides, so a tolerance would only hide a real difference.
agree(a::Integer, b::Integer) = a == b
agree(a::AbstractFloat, b::AbstractFloat) = isnan(a) ? isnan(b) : a === b

function main()
    write = "--write" in ARGS
    versions = setupR()
    rows = readfixture()
    println("Reference: ", versions)
    println("Fixture:   ", FIXTURE_PATH, " (", length(rows), " cases)\n")

    updated = similar(rows, 0)
    mismatches = String[]

    for row in rows
        reference = referencevalues(row.newick1, row.newick2)
        differing = [
            k for k in keys(reference) if !agree(reference[k], getproperty(row, k))
        ]

        if isempty(differing)
            println("  ok        ", row.case)
        else
            for k in differing
                push!(mismatches, "$(row.case): $k reference=$(repr(reference[k])) " *
                                  "committed=$(repr(getproperty(row, k)))")
            end
            println(write ? "  updated   " : "  MISMATCH  ", row.case,
                    " (", join(differing, ", "), ")")
        end
        push!(updated, merge(row, reference))
    end

    if write
        writefixture(updated)
        println("\nWrote ", FIXTURE_PATH)
        return 0
    end

    if isempty(mismatches)
        println("\nAll ", length(rows), " cases match the reference.")
        return 0
    end

    println("\n", length(mismatches), " value(s) differ from the reference:")
    foreach(m -> println("  ", m), mismatches)
    println("\nRe-run with --write to adopt the reference values, but read them first: " *
            "a disagreement means either this package or the fixture's stated derivation " *
            "is wrong.")
    return 1
end

exit(main())
