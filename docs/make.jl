using Nmrc
using Documenter
using Literate

DocMeta.setdocmeta!(Nmrc, :DocTestSetup, :(using Nmrc); recursive=true)

Literate.markdown("src/literate/flat.jl", "src/")

makedocs(;
    modules=[Nmrc],
    authors="Jiayi Zhang",
    sitename="Nmrc.jl",
    format=Documenter.HTML(;
        canonical="https://peakfind.github.io/Nmrc.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Tutorials" => "flat.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/peakfind/Nmrc.jl",
    devbranch="main",
)
