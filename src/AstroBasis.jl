module AstroBasis

using HDF5

# define an abstract Basis type; define a function to get the dimension for each basis type
abstract type AbstractAstroBasis end
dimension(basis::AbstractAstroBasis) = dimension(typeof(basis))


# split Basis types by dimension; specify the dimensionality for safety; correctly dispatch dimension for subtypes
abstract type RazorThinBasis <: AbstractAstroBasis end
dimension(::Type{<:RazorThinBasis}) = 2

abstract type SphericalBasis <: AbstractAstroBasis end
dimension(::Type{<:SphericalBasis}) = 3

# make dimension available
export dimension

# Documentation of the function needed in a basis
include("Utils/IO.jl")
export getparameters

# Documentation of the function needed in a basis
include("Basisdoc.jl")
export getUln,getDln,tabUl!,tabDl!

# bring in the Clutton-Brock (1973) spherical basis
include("CB73.jl")
export CB73Basis

# bring in the Hernquist & Ostriker (1992) disc basis
include("Hernquist.jl")
export HernquistBasis

# @IMPROVE, add Bessel basis

# @IMPROVE, add readers for EXP empirical bases

# bring in the Clutton-Brock (1972) disc basis
include("CB72.jl")
export CB72Basis

# bring in the Kalnajs (1976) disc basis
include("Kalnajs76.jl")
export K76Basis


end # module
