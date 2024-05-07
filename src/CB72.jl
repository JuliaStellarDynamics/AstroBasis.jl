"""
Radial basis elements from Clutton-Brock (1972)

@WARNING: Basis adapted to razor-thin discs, not spherical systems.
To use with azimuthal elements e^{i l phi} (not spherical harmonics).

@WARNING: Index starts at n=0. We therefore define nmax = nradial-1 the maximal
radial index (as lmax for azimuthal/harmonic index l).
"""

"""
    CB72Basis

Radial basis elements from Clutton-Brock (1972)
"""
struct CB72Basis <: AbstractAstroBasis

    name::String        # Basis name (default CB72)
    dimension::Int64     # Basis dimension (default 2)

    lmax::Int64         # Maximal harmonic/azimuthal index (starts at 0)
    nradial::Int64      # Number of radial basis elements (≥ 1)

    G::Float64          # Gravitational constant (default 1.0)
    rb::Float64         # Radial extension (default 1.0)

    tabPrefU::Array{Float64,2}  # Potential prefactors array
    tabPrefD::Array{Float64,2}  # Density prefactors array

    tabUl::Array{Float64,1}     # Potential elements value array
    tabDl::Array{Float64,1}     # Density elements value array

end

"""
    CB72Basis([name, dimension, lmax, nradial, G, rb])

creates a CB72Basis structure (and fill prefactors)
"""
function CB72Basis(;name::String="CB72", dimension::Int64=2,
                    lmax::Int64=0, nradial::Int64=0,
                    G::Float64=1., rb::Float64=1.)

    basis = CB72Basis(name,dimension,
                      lmax,nradial,
                      G,rb,
                      zeros(Float64,nradial,lmax+1),zeros(Float64,nradial,lmax+1), # Prefactors arrays
                      zeros(Int64,nradial),zeros(Float64,nradial)) # Elements value arrays

    FillPrefactors!(basis)

    return basis
end


##################################################
# Computation of the prefactors
##################################################
"""
    PrefactorsAzimuthalRecurrenceCB72(previous, l, n)

Gives the next adimensional prefactor (recurrence over azimuthal number l),
a^n_l = a^n_{l-1} * 2 * sqrt( 1 / [(n+2l)*(n+2l-1)] ),
given the previous one (a^n_{l-1}).

Initializing the recurrence at a^n_0 = sqrt(2), one gets the expected prefactor
a^n_l = 2^(l+1/2) * sqrt( n! / (n+2l)! ).
"""
function PrefactorsAzimuthalRecurrenceCB72(previous::Float64,
                                      l::Int64,
                                      n::Int64)::Float64
    return 2.0 * sqrt( 1.0 / ( (n+2.0*l)*(n+2.0*l-1.0) ) ) * previous
end

"""
    FillPrefactors!(basis::CB72Basis)

Clutton-Brock (1972) prefactors a^n_l = 2^(l+1/2) * sqrt( n! / (n+2l)! )
computed through the recurrence relation given by the function
prefactors_recurrence_l_CB72.

@IMPROVE precompute the a^n_l up to a given limit and save them in a file to read?
"""
function FillPrefactors!(basis::CB72Basis)

    lmax, nradial   = basis.lmax, basis.nradial
    G, rb           = basis.G, basis.rb
    tabPrefU, tabPrefD = basis.tabPrefU, basis.tabPrefD

    nmax = nradial - 1

    dimU = - sqrt(G / rb)                       # Potential basis element dimensional prefactor
    dimD = 1.0 / ( 2.0*pi*sqrt(G * (rb)^(3)) )  # Density basis element dimensional prefactor
    # Initialization
    for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        tabPrefU[n+1,1] = sqrt(2.0) * dimU
        tabPrefD[n+1,1] = sqrt(2.0) * dimD
    end
    # Recurrence
    for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        for l=1:lmax # Loop over the harmonic indices. ATTENTION, harmonic index starts at l=0 (here 0 has been initialized)
            tabPrefU[n+1,l+1] = PrefactorsAzimuthalRecurrenceCB72(tabPrefU[n+1,l],l,n)
            tabPrefD[n+1,l+1] = PrefactorsAzimuthalRecurrenceCB72(tabPrefD[n+1,l],l,n)
        end
    end

end


##################################################
# Computation of the potential basis elements
##################################################
"""
    ρCB72(x)

returns the parameter -1 <= ρ <= 1 for a given dimensionless radius x=r/rb.
"""
@inline function ρCB72(x::Float64)::Float64
    return (x*x - 1.0)/(x*x + 1.0) # Value of rho
end

"""
    UAzimuthalInitializationCB72(x, l)

Gives the (n=0, l) adimensional potential basis element value at a given position x = r/rb,
xi^0_l(x) = prod(2k - 1, k=1 -> l) / [  (1 + x^2)^{l + 1/2}   ].

@IMPROVE other function to give the initialization for all l simultaneously through recurrence?
"""
function UAzimuthalInitializationCB72(x::Float64,
                                              l::Int64)::Float64
    return prod(1:2:(2*l-1)) / (sqrt(1.0 + x*x) * (1.0 + x*x)^(l))
end

"""
    URadialReccurenceCB72(u^{n-2}, u^{n-1}, rho(x), l, n)

Gives the next adimensional potential basis element (recurrence over radial number n),
xi^n_l(x), given the 2 previous ones, xi^{n-1}_l(x) and xi^{n-2}_l(x).

Initializing the recurrence at xi^0_l given by the initialization from U_init_l_CB72,
one get the right next adimensional potential basis element from Clutton-Brock (1972).
"""
function URadialRecurrenceCB72(u0::Float64, u1::Float64,
                                ρ::Float64,
                                l::Int64, n::Int64)::Float64
    return ( 2.0 + ((2.0*l-1.0)/n) ) * ρ * u1 - ( 1.0 + ((2.0*l-1.0)/n) ) * u0
end

"""
    tabUl!(basis::CB72Basis, l, r[, forD])

For Clutton-Brock (1972) basis elements.
"""
function tabUl!(basis::CB72Basis,
                    l::Int64,r::Float64,
                    forD::Bool=false)

    nradial     = basis.nradial
    rb          = basis.rb
    tabPrefU    = basis.tabPrefU
    tabl        = forD ? basis.tabDl : basis.tabUl

    nmax = nradial - 1

    x   = r/rb        # Dimensionless radius
    xl  = x^(l)
    ρ = ρCB72(x)  # Value of the rescaled parameter rho

    #####
    # Recurrence on adimensional subpart (xi) of the potential basis elements
    #####
    # Initialization
    u0, u1 = 0.0, UAzimuthalInitializationCB72(x,l)  # n = -1, 0
    tabl[1] =  u1
    # Recurrence loop
    for n=1:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        u2 = URadialRecurrenceCB72(u0,u1,ρ,l,n)
        tabl[n+1] = u2
        u0, u1 = u1, u2
    end

    #####
    # Multiplying by the appropriate prefactors (except if nopref=true, then return only xi^n_l(x))
    #####
    if !(forD)
        for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
            tabl[n+1] *= tabPrefU[n+1,l+1] * xl
        end
    end
end

"""
    getUln(basis::CB72Basis, l, n, r[, forD])

For Clutton-Brock (1972) basis elements.
"""
function getUln(basis::CB72Basis,
                    l::Int64,n::Int64,r::Float64,
                    forD::Bool=false)

    rb          = basis.rb
    tabPrefU    = basis.tabPrefU

    x   = r/rb        # Dimensionless radius
    ρ   = ρCB72(x)  # Value of the rescaled parameter rho

    #####
    # Recurrence on adimensional subpart (xi) of the potential basis elements
    #####
    # Initialization
    uprev, uact = 0.0, UAzimuthalInitializationCB72(x,l)  # n = -1, 0
    # Recurrence loop
    for np=1:n # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        unext = URadialRecurrenceCB72(uprev,uact,ρ,l,np)
        uprev, uact = uact, unext
    end

    #####
    # Multiplying by the appropriate prefactor (except if nopref=true, then return only xi^n_l(x))
    #####
    if !(forD)
        uact *= tabPrefU[n+1,l+1] * x^(l)
    end

    return uact
end


##################################################
# Computation of the density basis elements
##################################################
"""
    tabDl!(basis::CB72Basis, l, r)

For Clutton-Brock (1972) basis elements, using tabUl!.
"""
function tabDl!(basis::CB72Basis,
                    l::Int64,r::Float64)
    #####
    # Compute potential basis elements at azimuthal number l+1 without prefactors
    # Stored in tabD
    #####
    tabUl!(basis,l+1,r,true)
    #####
    # Deduce density basis elements at azimuthal number l without prefactors
    #####
    nradial     = basis.nradial
    nmax        = nradial - 1
    rb          = basis.rb
    tabPrefD    = basis.tabPrefD
    tabDl       = basis.tabDl

    # For radial number n > 1.
    for n=nmax:-1:2 # Loop over the radial basis numbers. ATTENTION, index starts at n=0. In reverse order for no impact of the previous change.
        tabDl[n+1] -= tabDl[n-1]
    end

    #####
    # Adding prefactors
    #####
    x   = r/rb        # Dimensionless radius
    xl  = x^(l)
    for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        tabDl[n+1] *= tabPrefD[n+1,l+1] * xl
    end
end

"""
    getDln(basis::CB72Basis, l, n, r)

For Clutton-Brock (1972) basis elements, using getUln.
"""
function getDln(basis::CB72Basis,
                    l::Int64,n::Int64,r::Float64)
    #####
    # Compute potential basis elements n-2 and n at azimuthal number l+1 without prefactors
    #####
    u0, u2 = getUln(basis,l+1,n-2,r,true), getUln(basis,l+1,n,r,true)
    #####
    # # Deduce density basis elements at azimuthal number l without prefactors
    #####
    if n<2
        Dln = u2
    else
        Dln = u2 - u0
    end
    #####
    # Adding prefactors
    #####
    rb          = basis.rb
    tabPrefD    = basis.tabPrefD

    x   = r/rb        # Dimensionless radius
    xl  = x^(l)
    Dln *= tabPrefD[n+1,l+1] * xl

    return Dln
end
