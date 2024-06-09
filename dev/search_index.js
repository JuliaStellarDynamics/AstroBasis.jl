var documenterSearchIndex = {"docs":
[{"location":"common.html#Functions-common-to-all-bases","page":"Common Interface","title":"Functions common to all bases","text":"","category":"section"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"To seamlessly allow for switching between different bases, the bases are all organised around common functions. ","category":"page"},{"location":"common.html#Common-functions","page":"Common Interface","title":"Common functions","text":"","category":"section"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"The common functions are detailed in src/Basisdoc.jl. Any instantiated basis will result in a structure, with attributes name, dimension, lmax, nmax, rb (the basis scale), and G (the gravitational constant; unity by default). The following functions are also defined:","category":"page"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"getUln(basis,l,n,r) will calculate the basis potential at location r for harmonic order l and radial order n.","category":"page"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"getDln(basis,l,n,r) will calculate the basis density at location r for harmonic order l and radial order n.","category":"page"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"Overloaded function tabUl!(basis,l) will fill the structure attribute tabUl table with radial potential values from the basis up to nmax.","category":"page"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"Overloaded function tabDl!(basis,l) will fill the structure attribute tabUl table with radial density values from the basis up to nmax.","category":"page"},{"location":"common.html#Accessing-potential-functions","page":"Common Interface","title":"Accessing potential functions","text":"","category":"section"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"AstroBasis.getUln","category":"page"},{"location":"common.html#AstroBasis.getUln","page":"Common Interface","title":"AstroBasis.getUln","text":"getUln(basis,l,n,r)\n\nDefinition of the potential radial basis elements (l,n) at position r.\n\n@WARNING:: l=0 is index 1, n=0 is index 1.\n\n\n\n\n\ngetUln(basis::CB72Basis, l, n, r[, forD])\n\nFor Clutton-Brock (1972) basis elements.\n\n\n\n\n\ngetUln(basis::K76Basis, l, n, r[, forD])\n\nFor Kalnajs (1976) basis elements.\n\n\n\n\n\n","category":"function"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"AstroBasis.tabUl!","category":"page"},{"location":"common.html#AstroBasis.tabUl!","page":"Common Interface","title":"AstroBasis.tabUl!","text":"tabUl!(basis, l, r)\n\nCompute the basis potential elements for a given basis structure.\n\nArguments\n\nbasis: The structure containing basis parameters and prefactors. This should be a subtype of AbstractAstroBasis.\nl::Int64: The harmonic/azimuthal index.\nr::Float64: The radial coordinate.\n\nDescription\n\nThis function computes the basis elements for the provided structure. It updates the potential or density basis elements in-place within the given structure.\n\nImplementation Details\n\nThe function initializes the recurrence relation for the adimensional subpart of the basis elements.\nIt then iterates through the radial basis numbers, updating the basis elements using a recurrence relation.\nFinally, it multiplies the computed values by the appropriate prefactors.\n\nWarning\n\nl=0 is index 1, n=0 is index 1.\n\n\n\n\n\ntabUl!(basis::CB72Basis, l, r[, forD])\n\nCompute the Clutton-Brock (1972) basis elements for a given CB72Basis structure.\n\nArguments\n\nbasis::CB72Basis: The CB72Basis structure containing basis parameters and prefactors.\nl::Int64: The harmonic/azimuthal index.\nr::Float64: The radial coordinate.\nforD::Bool=false: A boolean flag indicating whether to use density prefactors. Default is false.\n\nDescription\n\nThis function computes the basis elements for the Clutton-Brock (1972) model. It updates the potential or density basis elements in-place within the provided CB72Basis structure.\n\nImplementation Details\n\nThe function initializes the recurrence relation for the adimensional subpart of the basis elements.\nIt then iterates through the radial basis numbers, updating the basis elements using a recurrence relation.\nFinally, it multiplies the computed values by the appropriate prefactors.\n\nExample\n\n# Assuming basis is an instance of CB72Basis with appropriate parameters\ntabUl!(basis, 2, 1.0)\n\nWarning\n\nEnsure that the tabPrefU and tabPrefD arrays in CB72Basis are correctly sized to avoid runtime errors.\n\n\n\n\n\ntabUl!(basis::K76Basis, l, r)\n\nFor Kalnajs (1976) basis elements.\n\n\n\n\n\n","category":"function"},{"location":"common.html#Accessing-density-functions","page":"Common Interface","title":"Accessing density functions","text":"","category":"section"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"AstroBasis.getDln","category":"page"},{"location":"common.html#AstroBasis.getDln","page":"Common Interface","title":"AstroBasis.getDln","text":"getDln(basis,l,n,r)\n\nDefinition of the density radial basis elements (l,n) at position r.\n\n@WARNING:: l=0 is index 1, n=0 is index 1.\n\n\n\n\n\ngetDln(basis::CB72Basis, l, n, r)\n\nFor Clutton-Brock (1972) basis elements, using getUln.\n\n\n\n\n\ngetDln(basis::K76Basis, l, n, r)\n\nFor Kalnajs (1976) basis elements.\n\n\n\n\n\n","category":"function"},{"location":"common.html","page":"Common Interface","title":"Common Interface","text":"AstroBasis.tabDl!","category":"page"},{"location":"common.html#AstroBasis.tabDl!","page":"Common Interface","title":"AstroBasis.tabDl!","text":"tabDl!(basis,l)\n\nCompute density table for a given l and r, and for 0 <= n <= nmax\n\nTo avoid repeated memory allocations, this overwrites the 'basis.tabDl' table\n\n@WARNING:: l=0 is index 1, n=0 is index 1.\n\n\n\n\n\ntabDl!(basis::CB72Basis, l, r)\n\nFor Clutton-Brock (1972) basis elements, using tabUl!.\n\n\n\n\n\ntabDl!(basis::K76Basis, l, r)\n\nFor Kalnajs (1976) basis elements.\n\n\n\n\n\n","category":"function"},{"location":"bases.html#Bases","page":"Bases","title":"Bases","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"","category":"page"},{"location":"bases.html#Basis-functions","page":"Bases","title":"Basis functions","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"AstroBasis currently offers four different types of bases.","category":"page"},{"location":"bases.html#3d-(Spherical)","page":"Bases","title":"3d (Spherical)","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"Clutton-Brock (1973), which is a match to the Plummer (1911) profile.\nHernquist & Ostriker (1992), which is a match to the Hernquist (1990) profile.","category":"page"},{"location":"bases.html#2d-(Razor-Thin-Discs)","page":"Bases","title":"2d (Razor-Thin Discs)","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"Clutton-Brock (1972), which is a 2d match to a Plummer profile.\nKalnajs (1976), which is a spiral basis.","category":"page"},{"location":"bases.html","page":"Bases","title":"Bases","text":"–","category":"page"},{"location":"bases.html#Spherical-bases","page":"Bases","title":"Spherical bases","text":"","category":"section"},{"location":"bases.html#Clutton-Brock-(1973)","page":"Bases","title":"Clutton-Brock (1973)","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"AstroBasis.CB73Basis","category":"page"},{"location":"bases.html#AstroBasis.CB73Basis","page":"Bases","title":"AstroBasis.CB73Basis","text":"CB73Basis <: SphericalBasis\n\nA structure for accessing the radial basis elements of Clutton-Brock (1973).\n\nFields\n\nname::String: Basis name (default CB73).\nlmax::Int64: Maximal harmonic/azimuthal index (starts at 0).\nnradial::Int64: Number of radial basis elements (≥ 1).\nG::Float64: Gravitational constant (default 1.0).\nrb::Float64: Radial extension (default 1.0).\ntabPrefU::Array{Float64,2}: Potential prefactors array.\ntabPrefD::Array{Float64,2}: Density prefactors array.\ntabUl::Array{Float64,1}: Potential elements value array.\ntabDl::Array{Float64,1}: Density elements value array.\n\nDescription\n\nThe CB73Basis structure is used to access the radial basis elements as defined by Clutton-Brock (1973). \n\nExample\n\n```julia using Random\n\nCreating an instance of CB73Basis\n\ncb73 = CB73Basis(     \"CB73\",      # name     10,          # lmax     5,           # nradial     1.0,         # G     1.0,         # rb     rand(10, 5), # tabPrefU     rand(10, 5), # tabPrefD     rand(5),     # tabUl     rand(5)      # tabDl )\n\nAccessing fields\n\nprintln(cb73.name) println(cb73.lmax) println(cb73.G)\n\nWarning\n\nEnsure that the prefactor arrays (tabPrefU and tabPrefD) and the element value arrays\n\n(tabUl and tabDl) have compatible dimensions to avoid runtime errors.\n\n\n\n\n\n","category":"type"},{"location":"bases.html#Hernquist-(1992)","page":"Bases","title":"Hernquist (1992)","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"AstroBasis.HernquistBasis","category":"page"},{"location":"bases.html#AstroBasis.HernquistBasis","page":"Bases","title":"AstroBasis.HernquistBasis","text":"HernquistBasis <: SphericalBasis\n\nA structure for interfacing with the radial basis elements from Hernquist & Ostriker (1992).\n\nFields\n\nname::String: Basis name (default Hernquist).\nlmax::Int64: Maximal harmonic/azimuthal index (starts at 0).\nnradial::Int64: Number of radial basis elements (≥ 1).\nG::Float64: Gravitational constant (default 1.0).\nrb::Float64: Radial extension (default 1.0).\ntabPrefU::Array{Float64,2}: Potential prefactors array.\ntabPrefD::Array{Float64,2}: Density prefactors array.\ntabUl::Array{Float64,1}: Potential elements value array.\ntabDl::Array{Float64,1}: Density elements value array.\n\nDescription\n\nThe HernquistBasis structure is used to access the radial basis elements as defined by Hernquist & Ostriker (1992).\n\n\n\n\n\n","category":"type"},{"location":"bases.html#Razor-thin-bases","page":"Bases","title":"Razor-thin bases","text":"","category":"section"},{"location":"bases.html#Clutton-Brock-(1972)","page":"Bases","title":"Clutton-Brock (1972)","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"AstroBasis.CB72Basis","category":"page"},{"location":"bases.html#AstroBasis.CB72Basis","page":"Bases","title":"AstroBasis.CB72Basis","text":"CB72Basis <: RazorThinBasis\n\nA structure for interfacing with the radial basis elements of Clutton-Brock (1972).\n\nFields\n\nname::String: Basis name (default CB72).\nlmax::Int64: Maximal harmonic/azimuthal index (starts at 0).\nnradial::Int64: Number of radial basis elements (≥ 1).\nG::Float64: Gravitational constant (default 1.0).\nrb::Float64: Radial extension (default 1.0).\ntabPrefU::Array{Float64,2}: Potential prefactors array.\ntabPrefD::Array{Float64,2}: Density prefactors array.\ntabUl::Array{Float64,1}: Potential elements value array.\ntabDl::Array{Float64,1}: Density elements value array.\n\nDescription\n\nThe CB72Basis structure is the interface to the radial basis elements as defined by Clutton-Brock (1972). \n\n\n\n\n\n","category":"type"},{"location":"bases.html#Kalnajs-(1976)","page":"Bases","title":"Kalnajs (1976)","text":"","category":"section"},{"location":"bases.html","page":"Bases","title":"Bases","text":"AstroBasis.K76Basis","category":"page"},{"location":"bases.html#AstroBasis.K76Basis","page":"Bases","title":"AstroBasis.K76Basis","text":"K76Basis <: RazorThinBasis\n\nA structure for interfacing with the radial basis elements of Kalnajs (1976).\n\nFields\n\nname::String: Basis name (default K76).\nlmax::Int64: Maximal harmonic/azimuthal index (starts at 0).\nnradial::Int64: Number of radial basis elements (≥ 1).\nG::Float64: Gravitational constant (default 1.0).\nrb::Float64: Radial extension (default 1.0).\nkKA::Int64: Basis index specific to Kalnajs (1976).\ntabPrefU::Array{Float64,2}: Potential prefactors array.\ntabPrefD::Array{Float64,2}: Density prefactors array.\ntabUl::Array{Float64,1}: Potential elements value array.\ntabDl::Array{Float64,1}: Density elements value array.\n\nDescription\n\nThe K76Basis structure is used to access the radial basis elements as defined by Kalnajs (1976). \n\nExample\n\nusing Random\n\n# Creating an instance of K76Basis\nk76 = K76Basis(\n    \"K76\",       # name\n    10,          # lmax\n    5,           # nradial\n    1.0,         # G\n    1.0,         # rb\n    2,           # kKA\n    rand(10, 5), # tabPrefU\n    rand(10, 5), # tabPrefD\n    rand(5),     # tabUl\n    rand(5)      # tabDl\n)\n\n# Accessing fields\nprintln(k76.name)\nprintln(k76.lmax)\nprintln(k76.kKA)\n\n\n\n\n\n","category":"type"},{"location":"quickstart.html#Quickstart-Example","page":"Quickstart","title":"Quickstart Example","text":"","category":"section"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"In this example code, we will make a figure of the radial basis elements from Clutton-Brock (1973).","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"import AstroBasis\nusing Plots\nusing LaTeXStrings","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"Now create the basis:","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"println(\"Creating the basis ... \")\nG, rb = 1., 1.\nltest, nradial = 2, 5\nbasis = AstroBasis.CB73Basis(lmax=ltest,nradial=nradial,G=G, rb=rb)","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"Define where to make the basis points:","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"# Points (rescaled radius)\nprintln(\"Compute basis values ... \")\nnx = 200\nrmin, rmax = 0., 3.\ntabx = collect(LinRange(rmin/basis.rb,rmax/basis.rb,nx))","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"Use the common function tabUl! to fill the table:","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"# Compute the values of the potential basis elements and store them\ntabU = zeros(nradial,nx) # Storage for the basis values\nfor j = 1:nx\n    # Compute all the basis elements values at a given location r (the result is stored in basis.tabUl)\n    AstroBasis.tabUl!(basis,ltest,tabx[j]*basis.rb)\n    # Store them in tabU\n    for i = 1:nradial\n        tabU[i,j] = basis.tabUl[i]\n    end\nend","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"And finally plot:","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"# Plot the curves\nprintln(\"Plotting ... \")\nlabels = reshape([\"n=\"*string(k) for k=1:nradial],(1,nradial)) #Need to be row\nplU=plot(tabx, transpose(tabU), title = \"Potential basis elements: Clutton-Brock (1973)\",label=labels)\nxlabel!(plU, L\"$r / r_{\\mathrm{b}}$\")\nylabel!(plU, L\"$U^{\\ell}_n (r)\\quad \\ell=$\"*string(ltest))\nsavefig(plU,\"CluttonBrock73.png\")\nprintln(\"The plot has been saved in the same folder as this example script under the name 'CluttonBrock73.png'.\")","category":"page"},{"location":"quickstart.html","page":"Quickstart","title":"Quickstart","text":"(Image: `Clutton-Brock (1973)`)","category":"page"},{"location":"index.html#AstroBasis.jl","page":"Home","title":"AstroBasis.jl","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"Bi-orthogonal bases for galactic dynamics in Julia","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"AstroBasis.jl implements several common basis sets (see 'Bases') in galactic dynamics, with a common interface (see 'Common Interfaces'). A quickstart is availabe, which will walk you through how to construct an image of the basis sets.","category":"page"},{"location":"index.html","page":"Home","title":"Home","text":"","category":"page"},{"location":"index.html#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"index.html","page":"Home","title":"Home","text":"For installation and first use test, follow the instructions on the github README.","category":"page"}]
}
