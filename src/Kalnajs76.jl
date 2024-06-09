"""
Radial basis elements from Kalnajs (1976)

@WARNING: Basis adapted to razor-thin discs, not spherical systems.
To use with azimuthal elements e^{i l phi} (not spherical harmonics).

@WARNING: Index starts at n=0. We therefore define nmax = nradial-1 the maximal
radial index (as lmax for azimuthal/harmonic index l).

@IMPROVE make rb only need to be set in one spot?

By default,
rb=1
G =1
"""

using SpecialFunctions
using Memoize

"""
K76Basis <: RazorThinBasis

A structure for interfacing with the radial basis elements of Kalnajs (1976).

# Fields
- `name::String`: Basis name (default K76).
- `lmax::Int64`: Maximal harmonic/azimuthal index (starts at 0).
- `nradial::Int64`: Number of radial basis elements (≥ 1).
- `G::Float64`: Gravitational constant (default 1.0).
- `rb::Float64`: Radial extension (default 1.0).
- `kKA::Int64`: Basis index specific to Kalnajs (1976).
- `tabPrefU::Array{Float64,2}`: Potential prefactors array.
- `tabPrefD::Array{Float64,2}`: Density prefactors array.
- `tabUl::Array{Float64,1}`: Potential elements value array.
- `tabDl::Array{Float64,1}`: Density elements value array.

# Description
The `K76Basis` structure is used to access the radial basis elements as defined
by Kalnajs (1976). 

# Example
```julia
using Random

# Creating an instance of K76Basis
k76 = K76Basis(
    "K76",       # name
    10,          # lmax
    5,           # nradial
    1.0,         # G
    1.0,         # rb
    2,           # kKA
    rand(10, 5), # tabPrefU
    rand(10, 5), # tabPrefD
    rand(5),     # tabUl
    rand(5)      # tabDl
)

# Accessing fields
println(k76.name)
println(k76.lmax)
println(k76.kKA)
```

"""
struct K76Basis <: RazorThinBasis

    name::String        # Basis name (default K76)

    lmax::Int64         # Maximal harmonic/azimuthal index (starts at 0)
    nradial::Int64      # Number of radial basis elements (≥ 1)

    G::Float64          # Gravitational constant (default 1.0)
    rb::Float64         # Radial extension (default 1.0)

    kKA::Int64          # Basis index

    tabPrefU::Array{Float64,2}  # Potential prefactors array
    tabPrefD::Array{Float64,2}  # Density prefactors array

    tabUl::Array{Float64,1}     # Potential elements value array
    tabDl::Array{Float64,1}     # Density elements value array

end

"""
    K76BasisCreate([name, dimension, lmax, nmax, G, rb])

creates a K76Basis structure
"""
function K76Basis(;name::String="K76",
                   lmax::Int64=0, nradial::Int64=1,
                   G::Float64=1., rb::Float64=1., kKA::Int64=1)

    basis = K76Basis(name,
                     lmax,nradial,
                     G,rb,kKA,
                     zeros(Float64,lmax+1,nradial),zeros(Float64,lmax+1,nradial), # Prefactors arrays
                     zeros(Int64,nradial),zeros(Float64,nradial)) # Elements value arrays

    FillPrefactors!(basis)

    return basis
end


##################################################
# Computation of the prefactors
##################################################
function PcoeffK76(k::Int64,l::Int64,n::Int64)
    return sqrt( (2*k+l+2*n+0.5) * gamma(2*k+l+n+0.5) * gamma(l+n+0.5) / (gamma(2*k+n+1) * (gamma(l+1))^(2) * gamma(n+1)) )
end
function ScoeffK76(k::Int64,l::Int64,n::Int64)
    return gamma(k + 1) / (pi * gamma(2*k+1) * gamma(k+0.5)) * 
            sqrt( (2*k+l+2*n+0.5) * gamma(2*k+n+1) * gamma(2*k+l+n+0.5) / (gamma(l+n+0.5) * gamma(n+1)) )
end
"""
    fill_prefactors!(basis::K76Basis)

Kalnajs (1976) prefactors
"""
function FillPrefactors!(basis::K76Basis)

    lmax, nradial  = basis.lmax, basis.nradial
    G, rb       = basis.G, basis.rb
    kKA         = basis.kKA
    tabPrefU, tabPrefD = basis.tabPrefU, basis.tabPrefD

    dimU = - sqrt(G / rb)                       # Potential basis element dimensional prefactor
    dimD = 1.0 / ( sqrt(G * (rb)^(3)) )  # Density basis element dimensional prefactor
    # Initialization
    for n=0:nradial-1 # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        for l=0:lmax # Loop over the harmonic indices. ATTENTION, harmonic index starts at l=0 
            tabPrefU[l+1,n+1] = dimU * PcoeffK76(kKA,l,n)
            tabPrefD[l+1,n+1] = (-1)^(n) * dimD * ScoeffK76(kKA,l,n)
        end
    end

end


##################################################
# Computation of the potential basis elements
##################################################
"""
    Pochhammer(x,n)

Pochhammer symbol for the rising factorial. From FastTransforms.jl package.
"""
@memoize function Pochhammer(x::Number,n::Integer)
    res = one(x)
    if n≥0
        for i=0:n-1
            res *= x+i
        end
    else
        res /= pochhammer(x+n,-n)
    end
    return res
end

@memoize function alphaK76(k::Int64, l::Int64, n::Int64,
                    i::Int64, j::Int64)

    return Pochhammer(-k,i) * Pochhammer(l+0.5,i) * Pochhammer(2*k+l+n+0.5,j) * Pochhammer(i+l+0.5,j) * Pochhammer(-n,j) /
            ( Pochhammer(l+1,i) * Pochhammer(1,i) * Pochhammer(l+i+1,j) * Pochhammer(l+0.5,j) * Pochhammer(1,j) )
end

"""
    tabUl!(basis::K76Basis, l, r)

For Kalnajs (1976) basis elements.
"""
function tabUl!(basis::K76Basis,
                    l::Int64,r::Float64)

    nradial     = basis.nradial
    nmax        = nradial - 1
    rb          = basis.rb
    kKA         = basis.kKA
    tabPrefU    = basis.tabPrefU
    tabUl       = basis.tabUl

    x   = r/rb        # Dimensionless radius

    fill!(tabUl,0.0)

    if x > 1.0
    else
        for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
            for i=0:kKA
                for j=0:n
                    tabUl[n+1] += alphaK76(kKA,l,n,i,j) * (x)^(2*i+2*j)
                end
            end
        end
        for n=0:nmax
            tabUl[n+1] *= tabPrefU[l+1,n+1] * (x)^(l)
        end
    end

end

"""
    getUln(basis::K76Basis, l, n, r[, forD])

For Kalnajs (1976) basis elements.
"""
function getUln(basis::K76Basis,
                    l::Int64,n::Int64,r::Float64)

    rb          = basis.rb
    kKA         = basis.kKA
    tabPrefU    = basis.tabPrefU

    x   = r/rb        # Dimensionless radius

    res = 0.0
    if x > 1.0
    else
        for i=0:kKA
            for j=0:n
                res += alphaK76(kKA,l,n,i,j) * (x)^(2*i+2*j)
            end
        end
        res *= tabPrefU[l+1,n+1] * (x)^(l)
    end

    return res
end


##################################################
# Computation of the density basis elements
##################################################
@memoize function betaK76(k::Int64, l::Int64, n::Int64,
                    j::Int64)

    return Pochhammer(2*k+l+n+0.5,j) * Pochhammer(k+1,j) * Pochhammer(-n,j) / 
                ( Pochhammer(2*k+1,j) * Pochhammer(k+0.5,j) * Pochhammer(1,j) )
end

"""
    tabDl!(basis::K76Basis, l, r)

For Kalnajs (1976) basis elements.
"""
function tabDl!(basis::K76Basis,
                    l::Int64,r::Float64)

    nradial     = basis.nradial
    nmax        = nradial - 1
    rb          = basis.rb
    kKA         = basis.kKA
    tabPrefD    = basis.tabPrefD
    tabDl       = basis.tabDl

    x   = r/rb        # Dimensionless radius
    
    fill!(tabDl,0.0)

    if x <= 1.0
        for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
            for j=0:n
                tabDl[n+1] += betaK76(kKA,l,n,j) * (1.0-x*x)^(j)
            end
        end
        for n=0:nmax
            tabDl[n+1] *= tabPrefD[l+1,n+1] * (1.0-x*x)^(kKA-0.5) * (x)^(l)
        end
    end
end

"""
    getDln(basis::K76Basis, l, n, r)

For Kalnajs (1976) basis elements.
"""
function getDln(basis::K76Basis,
                    l::Int64,n::Int64,r::Float64)
    
    rb          = basis.rb
    kKA         = basis.kKA
    tabPrefD    = basis.tabPrefD

    x   = r/rb        # Dimensionless radius

    res = 0.0
    if x <= 1.0
        for j=0:n
            res += betaK76(kKA,l,n,j) * (1.0-x*x)^(j)
        end
        res *= tabPrefD[l+1,n+1] * (1.0-x*x)^(kKA-0.5) * (x)^(l)
    end

    return res
end
