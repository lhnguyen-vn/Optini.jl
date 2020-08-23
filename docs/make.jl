# Remove when package is registered
push!(LOAD_PATH, "../src")

using Optini
using Documenter

makedocs(;
    modules=[Optini],
    authors="Long Nguyen",
    repo="https://github.com/lhnguyen-vn/Optini.jl/blob/{commit}{path}#L{line}",
    sitename="Optini.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://lhnguyen-vn.github.io/Optini.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/lhnguyen-vn/Optini.jl.git",
    devbranch = "main"
)
