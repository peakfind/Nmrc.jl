using Nmrc
using Documenter
using Literate

DocMeta.setdocmeta!(Nmrc, :DocTestSetup, :(using Nmrc); recursive=true)

source_path = joinpath(@__DIR__, "src", "literate", "flat.jl")
output_path = joinpath(@__DIR__, "src")
Literate.markdown(source_path, output_path)

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
