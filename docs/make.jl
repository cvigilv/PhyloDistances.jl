using DemoCards
using Documenter
using DocumenterCitations
using DocumenterCodeBlocks
using DocumenterVitepress
using PhyloDistances

# DocumenterVitepress publishes only the asset files in `src/assets` whose names contain
# "logo", plus `favicon.ico`. A file named exactly `logo.svg` or `logo.png` would also be
# placed in the nav bar, which the site leaves without a logo. The logo set itself lives in
# `assets/logo`.
logodir = joinpath(@__DIR__, "assets", "logo")
siteassets = mkpath(joinpath(@__DIR__, "src", "assets"))
for (src, dst) in [
        "phylodistances-lockup.svg" => "logo-lockup.svg",
        "phylodistances-lockup-dark.svg" => "logo-lockup-dark.svg",
        "favicon.ico" => "favicon.ico",
    ]
    cp(joinpath(logodir, src), joinpath(siteassets, dst); force = true)
end

examples, finish_examples, _ = makedemos("demos")
bibliography = CitationBibliography(joinpath(@__DIR__, "src", "references.bib"))

makedocs(;
    modules = [PhyloDistances],
    authors = "Carlos Vigil-Vásquez and contributors",
    sitename = "PhyloDistances.jl",
    repo = "https://github.com/cvigilv/PhyloDistances.jl",
    format = MarkdownVitepress(;
        repo = "https://github.com/cvigilv/PhyloDistances.jl",
        devbranch = "main",
        devurl = "dev",
    ),
    pages = [
        "Home" => "index.md",
        "API reference" => "api.md",
        examples,
    ],
    plugins = [
        CodeBlocks(),
        bibliography,
    ],
    checkdocs = :exports,
)

finish_examples()

DocumenterVitepress.deploydocs(;
    repo = "github.com/cvigilv/PhyloDistances.jl.git",
    devbranch = "main",
    push_preview = true,
)
