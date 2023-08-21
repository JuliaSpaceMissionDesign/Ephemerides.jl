using Documenter, Ephemerides

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="Ephemerides.jl",
    modules=[Ephemerides],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",
        "API" => [
            "API" => "api.md"
        ],
    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/Ephemerides.jl", branch="gh-pages"
)
