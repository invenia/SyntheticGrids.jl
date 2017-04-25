using Documenter, SyntheticGrids

makedocs(
    modules = [SyntheticGrids],
    format = :html,
    pages = [
        "Home" => "index.md",
        "Types" => "Types.md",
        "Model" => "Model.md",
        "Functions" => [
		"Functions.md",
		"Private.md",
	],
    ],
    repo = "https://github.com/invenia/SyntheticGrids.jl/blob/{commit}{path}#L{line}",
    sitename = "SyntheticGrids.jl",
    authors = "Invenia Technical Computing",
    assets = ["assets/invenia.css"],
)

deploydocs(
    repo = "github.com/invenia/SyntheticGrids.jl",
    deps   = Deps.pip("python-markdown-math"),
    julia  = "0.6",
    target = "build",
    make = nothing,
)

