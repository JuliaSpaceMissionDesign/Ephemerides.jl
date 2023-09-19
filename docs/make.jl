using Documenter, Ephemerides

const CI = get(ENV, "CI", "false") == "true"

makedocs(;
    authors="Julia Space Mission Design Development Team",
    sitename="Ephemerides.jl",
    modules=[Ephemerides],
    format=Documenter.HTML(; prettyurls=CI, highlights=["yaml"], ansicolor=true),
    pages=[
        "Home" => "index.md",
        
        "Tutorials" => [
            "Loading Kernels" => "tutorials/load.md",
            "Reading Ephemeris Data" => "tutorials/position.md",
            "Kernels Inspection" => "tutorials/inspect.md"
        ],

        "Benchmarks" => "benchmarks.md",

        "API" => [
            "Public API" => "api/api.md",
            "Low-level API" => "api/lapi.md"
        ]
    ],
    clean=true,
)

deploydocs(;
    repo="github.com/JuliaSpaceMissionDesign/Ephemerides.jl", branch="gh-pages"
)
