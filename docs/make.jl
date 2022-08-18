using Documenter, MimiBRICK

makedocs(
	doctest = false,
    modules = [MimiBRICK],
	sitename = "MimiBRICK.jl",
	pages = [
		"Home" => "index.md",
		"Installation" => "src/installation.md",
		"Calibratoin" => "src/calibration.md",
		"Examples" => "src/examples.md",
	],
	format = Documenter.HTML(prettyurls = get(ENV, "JULIA_NO_LOCAL_PRETTY_URLS", nothing) === nothing)
)

deploydocs(
    repo = "github.com/raddleverse/MimiBRICK.jl.git",
)
