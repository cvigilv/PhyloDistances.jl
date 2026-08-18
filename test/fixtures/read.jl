# Reading and writing `rf_quartet.tsv`.
#
# Shared by the portable test suite, which checks the committed values, and by
# `validation/fixture.jl`, which regenerates them from the R references. One reader means
# the file's format is defined in exactly one place.

const FIXTURE_PATH = joinpath(@__DIR__, "rf_quartet.tsv")

const FIXTURE_COLUMNS = [
    "case", "newick1", "newick2",
    "rf", "rf_normalized", "quartet", "quartet_normalized",
    "provenance",
]

const FIXTURE_HEADER = """
# Expected Robinson-Foulds and quartet distances for named tree pairs.
#
# The numeric columns are generated from TreeDist and Quartet by
# `validation/fixture.jl`, which also re-checks them against those packages; the
# `provenance` column gives the independent derivation that corroborates each row. See
# README.md in this directory for the column meanings.
#
# Tab-separated. Lines beginning with "#" and blank lines are comments; the first
# remaining line names the columns.
"""

"""
    readfixture(path = FIXTURE_PATH) -> Vector{<:NamedTuple}

One named tuple per case, with the two trees left as Newick strings.

Throws if the column names are not the expected ones or a row has the wrong field count,
so a malformed file fails immediately rather than yielding rows with shifted values.
"""
function readfixture(path = FIXTURE_PATH)
    rows = NamedTuple[]
    columns = nothing
    for line in eachline(path)
        (isempty(strip(line)) || startswith(line, '#')) && continue
        fields = split(line, '\t')

        if columns === nothing
            columns = fields
            columns == FIXTURE_COLUMNS || error(
                "$path: columns are $columns, expected $FIXTURE_COLUMNS"
            )
            continue
        end

        length(fields) == length(columns) || error(
            "$path: case $(first(fields)) has $(length(fields)) fields, " *
            "expected $(length(columns))"
        )
        push!(rows, (
            case = String(fields[1]),
            # `readnw` rejects a `SubString`, reporting a missing semicolon.
            newick1 = String(fields[2]),
            newick2 = String(fields[3]),
            rf = parse(Int, fields[4]),
            rf_normalized = parse(Float64, fields[5]),
            quartet = parse(Int, fields[6]),
            quartet_normalized = parse(Float64, fields[7]),
            provenance = String(fields[8]),
        ))
    end
    columns === nothing && error("$path: no header line")
    return rows
end

"""
    writefixture(rows, path = FIXTURE_PATH)

Render `rows` back to the fixture file.

Numbers are written with Julia's shortest round-tripping representation, so reading the
file returns the very `Float64` that was written.
"""
function writefixture(rows, path = FIXTURE_PATH)
    open(path, "w") do io
        print(io, FIXTURE_HEADER)
        println(io, join(FIXTURE_COLUMNS, '\t'))
        for row in rows
            any(contains('\t'), (row.case, row.newick1, row.newick2, row.provenance)) &&
                error("$(row.case): a tab in a field would corrupt the file")
            println(io, join((
                row.case, row.newick1, row.newick2,
                row.rf, row.rf_normalized, row.quartet, row.quartet_normalized,
                row.provenance,
            ), '\t'))
        end
    end
    return path
end
