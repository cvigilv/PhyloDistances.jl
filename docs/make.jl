using DemoCards
using Documenter
using DocumenterCitations
using DocumenterCodeBlocks
using DocumenterVitepress
using PhyloDistances

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
