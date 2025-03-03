using Nmrc
using Documenter

DocMeta.setdocmeta!(Nmrc, :DocTestSetup, :(using Nmrc); recursive=true)

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
    ],
)

deploydocs(;
    repo="github.com/peakfind/Nmrc.jl",
    devbranch="main",
)
