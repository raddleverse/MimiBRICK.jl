using MimiBRICK
using Documenter
using Missings,DataFrames

makedocs(
	doctest = false,
    modules = [MimiBRICK],
	sitename = "MimiBRICK.jl Documentation",
	pages = ["Home" => "index.md",
	         "Installation and Examples" => "install_and_examples.md"
	],
	format = Documenter.HTML(prettyurls = get(ENV, "JULIA_NO_LOCAL_PRETTY_URLS", nothing) === nothing)
)

deploydocs(
    repo = "github.com/raddleverse/MimiBRICK.jl.git",
)
