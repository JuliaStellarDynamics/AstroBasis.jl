module AstroBasis

using HDF5

# Documentation of the function needed in a basis
include("Basisdoc.jl")

# bring in the Clutton-Brock (1973) spherical basis
include("CB73.jl")

# @IMPROVE, add Bessel basis

# bring in the Hernquist & Ostriker (1992) disc basis
include("Hernquist.jl")

# bring in the Clutton-Brock (1972) disc basis
include("CB72.jl")

# bring in the Kalnajs (1976) disc basis
include("Kalnajs76.jl")

# @IMPROVE, add readers for EXP empirical bases


# make a generic Basis data type
BasisType = Union{CB73Basis,HernquistBasis,CB72Basis,K76Basis}


end # module
