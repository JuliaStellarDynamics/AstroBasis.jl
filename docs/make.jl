# Adding the package src to the load path
push!(LOAD_PATH,"../src/")

using Documenter, AstroBasis

makedocs(sitename = "AstroBasis.jl",
         pages=[
                "Home" => "index.md",
                "Quickstart" => "quickstart.md",
                "Common Interface" => "common.md",
                "Bases" => "bases.md"
               ],
         format = Documenter.HTML(prettyurls=false))

deploydocs(repo="github.com/JuliaStellarDynamics/AstroBasis.jl",devbranch="documentation")