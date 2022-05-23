module AstroBasis

using HDF5

# Documentation of the function needed in a basis
include("Basisdoc.jl")

# bring in the Clutton-Brock basis
include("CB73.jl")

# @IMPROVE, add Bessel basis

# @IMPROVE, add Hernquist basis

# @IMPROVE, add Clutton-Brock (1972) disc basis
include("CB72.jl")

# @IMPROVE, add Kalnajs (1976) disc basis

# @IMPROVE, add readers for EXP empirical bases

end # module
