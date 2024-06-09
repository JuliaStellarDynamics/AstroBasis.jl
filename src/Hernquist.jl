"""
Radial basis elements from Hernquist & Ostriker (1992)

"""

# set the data path for the basis prefactor elements
data_path_hernquist() = abspath(joinpath(@__DIR__, "tables", "data_HO92_lmax_50_nmax_200.h5"))


"""
HernquistBasis <: SphericalBasis

A structure for interfacing with the radial basis elements from Hernquist & Ostriker (1992).

# Fields
- `name::String`: Basis name (default Hernquist).
- `lmax::Int64`: Maximal harmonic/azimuthal index (starts at 0).
- `nradial::Int64`: Number of radial basis elements (≥ 1).
- `G::Float64`: Gravitational constant (default 1.0).
- `rb::Float64`: Radial extension (default 1.0).
- `tabPrefU::Array{Float64,2}`: Potential prefactors array.
- `tabPrefD::Array{Float64,2}`: Density prefactors array.
- `tabUl::Array{Float64,1}`: Potential elements value array.
- `tabDl::Array{Float64,1}`: Density elements value array.

# Description
The `HernquistBasis` structure is used to access the radial basis elements as defined
by Hernquist & Ostriker (1992).

"""
struct HernquistBasis <: SphericalBasis

    name::String         # Basis name (default Hernquist)

    lmax::Int64         # Maximal harmonic/azimuthal index (starts at 0)
    nradial::Int64      # Number of radial basis elements (≥ 1)

    G::Float64       # Gravitational constant (default 1.0)
    rb::Float64      # Radial extension (default 1.0)

    tabPrefU::Array{Float64,2}  # Potential prefactors array
    tabPrefD::Array{Float64,2}  # Density prefactors array

    tabUl::Array{Float64,1}     # Potential elements value array
    tabDl::Array{Float64,1}     # Density elements value array

end


"""
HernquistBasis([name, dimension, lmax, nmax, G, rb, filename])

Create a HernquistBasis structure and fill prefactors.

# Arguments
- `name::String="Hernquist"`: Basis name.
- `lmax::Int64=0`: Maximal harmonic/azimuthal index (starts at 0).
- `nradial::Int64=1`: Number of radial basis elements (≥ 1).
- `G::Float64=1.0`: Gravitational constant.
- `rb::Float64=1.0`: Radial extension.
- `filename::String=data_path_hernquist()`: File path to Hernquist prefactors data.

# Description
The `HernquistBasis` constructor function creates a new instance of the `HernquistBasis` structure.
It initializes the basis with provided parameters and fills the potential and density prefactor arrays
by reading from a file specified by `filename`.

# Example
```julia
# Creating a HernquistBasis with default parameters
hernquist_basis = HernquistBasis()

# Creating a HernquistBasis with custom parameters
custom_basis = HernquistBasis(name="CustomBasis", lmax=10, nradial=5, G=2.0, rb=2.0)
```

"""
function HernquistBasis(;name::String="Hernquist",
                         lmax::Int64=0, nradial::Int64=1,
                         G::Float64=1., rb::Float64=1.,
                         filename::String=data_path_hernquist())

    tabPrefU, tabPrefD = ReadFillHernquistPrefactors(lmax,nradial,rb,G,filename)
    basis = HernquistBasis(name,
                           lmax,nradial,
                           G,rb,
                           tabPrefU,tabPrefD, # Prefactors arrays
                           zeros(Int64,nradial),zeros(Float64,nradial)) # Elements value arrays

    return basis
end


"""ReadFillHernquistPrefactors(lmax,nradial[,rb,G,precomputed_filename])

reads a table of pre-computed prefactors for Gegenbauer functions
comes loaded with a pre-computed large table of prefactors (probably more than you need!)

"""
function ReadFillHernquistPrefactors(lmax::Int64,nradial::Int64,
                                     rb::Float64=1.,G::Float64=1.,
                                     filename::String=data_path_hernquist())

    tabαHernquist = h5read(filename,"tab_alphalnp") # Reading the prefactors α_lnp
    tabβHernquist  = h5read(filename,"tab_betalnp")  # Reading the prefactors beta_lnp


    # generate blank tables: set up as constants; may only be computed once per session
    tabPrefHernquist_Ulnp = zeros(Float64,lmax+1,nradial) # Table of the prefactors of the POTENTIAL basis functions, i.e. the Ulnp !! ATTENTION, size is lmax+1 as l=0 exists
    tabPrefHernquist_Dlnp = zeros(Float64,lmax+1,nradial) # Table of the prefactors of the DENSITY   basis functions, i.e. the Dlnp !! ATTENTION, size is lmax+1 as l=0 exists

    for l=0:lmax # Loop over the harmonic indices. ATTENTION, harmonic index starts at l=0
        for n=0:nradial-1 # Loop over the radial basis numbers
            α = tabαHernquist[l+1,n+1]  # Reading the value of α. l,np starts at 0
            β  = tabβHernquist[l+1,n+1]   # Reading the value of β.  l,np starts at 0

            A = sqrt(G/rb)*α            # Value of the prefactor A_n^\ell
            B = 1.0/(sqrt(G)*rb^(5/2))*β # Value of the prefactor B_n^\ell

            tabPrefHernquist_Ulnp[l+1,n+1] = A  # Filling in the array. l,np starts at 0
            tabPrefHernquist_Dlnp[l+1,n+1] = B  # Filling in the array. l,np starts at 0
        end
    end

    return tabPrefHernquist_Ulnp, tabPrefHernquist_Dlnp
end

"""
Function that returns the parameter -1 <= xi <= 1 for a given dimensionless radius x=r/Rbasis
"""
function rhoHernquist(x::Float64)
    return (x - 1.0)/(x + 1.0) # Value of xi
end



"""ClnHernquist(α,n,rho)
Definition of the Gegenbauer polynomials
These coefficients are computed through an upward recurrence
@IMPROVE compute all the basis elements (n)_{1<=n<=nmax} at once
"""
function ClnHernquist(α::Float64,n::Int64,xi::Float64)

    v0 = 1.0 # Initial value for n=0
    if (n == 0)
        return v0 # No need for a recurrence for n=0
    end

    v1 = 2.0*α*xi # Initial value for n=1
    if (n == 1)
        return v1 # No need for a recurrence for n=1
    end

    ic = 2 # Iteration counter that gives the index of the value that is about to be computed
    v = 0.0 # Initialisation of the temporary variable

    while (ic <= n) # Applying the recurrence as many times as needed
        v = (2.0*(ic+α-1.0)*xi*v1 - (ic+2.0*α-2.0)*v0)/(ic) # Applying the recurrence
        v0, v1 = v1, v # Updating the temporary variables
        ic += 1 # Updating the counter of iteration
    end

    return v # Output of the value
end

"""UlnpHernquist(lharmonic,np,r,prefactor_table[,rb])
Definition of the potential radial basis elements from Clutton-Brock (1973)

Be careful: l=0 is index 1, np=0 is index 1.
"""
function UlnpHernquist(l::Int64,np::Int64,r::Float64,tabPrefHernquist_Ulnp::Matrix{Float64},rb::Float64=1.)
    pref = tabPrefHernquist_Ulnp[l+1,np+1]              # Value of the prefactor. l,np start l=0
    x    = r/rb                                    # Dimensionless radius
    rho  = rhoHernquist(x)                              # Value of the rescaled parameter rho
    valR = ((x/(1.0+x^(2)))^(l))/(sqrt(1.0+x^(2))) # Value of the multipole factor
    valC = ClnHernquist(l+1.0,np,rho)                 # Value of the Gegenbauer polynomials
    res  = pref*valR*valC                          # Value of the radial function
    return res
end

function getUln(basis::HernquistBasis,l::Int64,np::Int64,r::Float64)
    return UlnpHernquist(l,np,r,basis.tabPrefU,basis.rb)
end


"""DlnpHernquist(lharmonic,np,r,prefactor_table[,rb])
Definition of the density radial basis elements from Clutton-Brock (1973)

Be careful: l=0 is index 1, np=0 is index 1.
"""
function DlnpHernquist(l::Int64,np::Int64,r::Float64,tabPrefHernquist_Dlnp::Matrix{Float64},rb::Float64=1.)
    pref = tabPrefHernquist_Dlnp[l+1,np+1]                # Value of the prefactor. ATTENTION, l starts at l = 0
    x    = r/rb                                      # Dimensionless radius
    rho  = rhoHernquist(x)                                # Value of the rescaled parameter rho
    valR = ((x/(1.0+x^(2)))^(l))/((1.0+x^(2))^(5/2)) # Value of the multipole factor
    valC = ClnHernquist(l+1.0,np,rho)                   # Value of the Gegenbauer polynomials
    res  = pref*valR*valC
    return res
end


function getDln(basis::HernquistBasis,l::Int64,np::Int64,r::Float64)
    return DlnpHernquist(l,np,r,basis.tabPrefD,basis.rb)
end


"""tabUlnpHernquist!(lharmonic,radius,tabUlnp,nmax,prefactor_table[,rb])

Compute Hernquist potential table for a given l and r, and for 1 <= np <= nmax, nearly simultaneously

To avoid repeated memory allocations, this overwrites a table already in the struct.

@IMPROVE, create outs for nmax<3?
"""
function tabUlnpHernquist!(l::Int64,r::Float64,
                      tabUlnp::Array{Float64,1},
                      nmax::Int64,
                      tabPrefHernquist_Ulnp::Matrix{Float64},
                      rb::Float64=1.)

    x     = r/rb                                    # Dimensionless radius
    rho   = rhoHernquist(x)                              # Value of the rescaled parameter rho
    valR  = ((x/(1.0+x^(2)))^(l))/(sqrt(1.0+x^(2))) # Value of the multipole factor

    α = l+1.0 # Value of α, the index of the Gegenbauer polynomials

    # set initial values for recurrence
    v0 = 1.0           # Initial value of the Gegenbauer polynomials for n=0
    tabUlnp[1] = tabPrefHernquist_Ulnp[l+1,1]*valR*v0 # Filling in the value for np=1. ATTENTION, l starts at l=0

    v1 = 2.0*α*rho # Initial value of the Gegenbauer polynomials for n=1
    tabUlnp[2] = tabPrefHernquist_Ulnp[l+1,2]*valR*v1 # Filling in the value for np=2. ATTENTION, l starts at l=0

    for np=2:nmax-1 # Loop over remaining the radial indices
        v = (2.0*(np+α-1.0)*rho*v1 - (np+2.0*α-2.0)*v0)/(np) # Applying the recurrence
        v0, v1 = v1, v # Updating the temporary variables
        tabUlnp[np+1] = tabPrefHernquist_Ulnp[l+1,np+1]*valR*v # Filling in the value for np. ATTENTION, l starts at l=0
    end
end


function tabUl!(basis::HernquistBasis,l::Int64,r::Float64)

    tabUlnpHernquist!(l,r,basis.tabUl,basis.nradial,basis.tabPrefU,basis.rb)
    #return UlnpHernquist(l,np,r,basis.tabPrefU,basis.rb)
end


"""tabDlnpHernquist!(lharmonic,radius,tabD,nmax,prefactor_table[,rb])

Compute Hernquist density table for a given l and r, and for 1 <= np <= nmax, nearly simultaneously

To avoid repeated memory allocations, this overwrites a table already in the struct.

@IMPROVE, create outs for nmax<3?
"""
function tabDlnpHernquist!(l::Int64,r::Float64,
                      tabDlnp::Array{Float64,1},
                      nmax::Int64,
                      tabPrefHernquist_Dlnp::Matrix{Float64},
                      rb::Float64=1.)

    x     = r/rb                                      # Dimensionless radius
    rho   = rhoHernquist(x)                                # Value of the rescaled parameter rho
    valR  = ((x/(1.0+x^(2)))^(l))/((1.0+x^(2))^(5/2)) # Value of the multipole factor

    α = l+1.0 # Value of α, the index of the Gegenbauer polynomials

    # set initial values for recurrence
    v0 = 1.0           # Initial value of the Gegenbauer polynomials for n=0
    tabDlnp[1] = tabPrefHernquist_Dlnp[l+1,1]*valR*v0 # Filling in the value for np=1. ATTENTION, l starts at l=0

    v1 = 2.0*α*rho # Initial value of the Gegenbauer polynomials for n=1
    tabDlnp[2] = tabPrefHernquist_Dlnp[l+1,2]*valR*v1 # Filling in the value for np=2. ATTENTION, l starts at l=0

    for np=2:nmax-1 # Loop over remaining the radial indices
        v = (2.0*(np+α-1.0)*rho*v1 - (np+2.0*α-1.0)*v0)/(np) # Applying the recurrence
        v0, v1 = v1, v # Updating the temporary variables
        tabDlnp[np+1] = tabPrefHernquist_Dlnp[l+1,np+1]*valR*v # Filling in the value for np. ATTENTION, l starts at l=0
    end
end


function tabDl!(basis::HernquistBasis,l::Int64,r::Float64)

    tabDlnpHernquist!(l,r,basis.tabDl,basis.nradial,basis.tabPrefD,basis.rb)
end
