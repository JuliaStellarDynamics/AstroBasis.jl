"""
Radial basis elements from Clutton-Brock (1972)

@WARNING: Basis adapted to razor-thin discs, not spherical systems.
To use with azimuthal elements e^{i l phi} (not spherical harmonics).

@WARNING: Index starts at n=0. We therefore define nmax = nradial-1 the maximal 
radial index (as lmax for azimuthal/harmonic index l).

@IMPROVE make rb only need to be set in one spot?

By default,
rb=1
G =1
"""

"""
    structCB72Basis_type

Radial basis elements from Clutton-Brock (1972)
"""
struct structCB72Basis_type

    name::String        # Basis name (default CB72)
    dimension::Int64     # Basis dimension (default 2)

    lmax::Int64         # Maximal harmonic/azimuthal index (starts at 0)
    nmax::Int64         # Maximal radial index (starts at 0)

    G::Float64          # Gravitational constant (default 1.0)
    rb::Float64         # Radial extension (default 1.0)

    tabPrefU::Array{Float64,2}  # Potential prefactors array
    tabPrefD::Array{Float64,2}  # Density prefactors array

    tabUl::Array{Float64,1}     # Potential elements value array
    tabDl::Array{Float64,1}     # Density elements value array

end

"""
    CB72Basis_create([name, dimension, lmax, nmax, G, rb])

Create a structCB72Basis_type structure

By default, 
name="CB72", dimension=2,
lmax=0, nmax=0,
G=1., rb=1.
"""
function CB72Basis_create(name::String="CB72", dimension::Int64=2,
                            lmax::Int64=0, nmax::Int64=0,
                            G::Float64=1., rb::Float64=1.)
    return structMultipole_type(name,dimension, 
                                lmax,nmax,
                                G,rb,
                                zeros(Float64,lmax+1,nmax+1),zeros(Float64,lmax+1,nmax+1), # Prefactors arrays
                                zeros(Int64,nmax+1),zeros(Float64,nmax+1)) # Elements value arrays
end


##################################################
# Computation of the prefactors 
##################################################
"""
    prefactors_recurrence_l_CB72(previous, l, n)

Gives the next adimensional prefactor (recurrence over azimuthal number l),
a^n_l = a^n_{l-1} * 2 * sqrt( 1 / [(n+2l)*(n+2l-1)] ),
given the previous one (a^n_{l-1}).

Initializing the recurrence at a^n_0 = sqrt(2), one gets the expected prefactor
a^n_l = 2^(l+1/2) * sqrt( n! / (n+2l)! ).
"""
function prefactors_recurrence_l_CB72(previous::Float64,
                                l::Int64,
                                n::Int64)
    return 2.0 * sqrt( 1.0 / ( (n+2.0*l)*(n+2.0*l-1.0) ) ) * previous
end

"""
    fill_prefactors!(basis::structCB72Basis_type)

Clutton-Brock (1972) prefactors a^n_l = 2^(l+1/2) * sqrt( n! / (n+2l)! )
computed through the recurrence relation given by the function
prefactors_recurrence_l_CB72.

@IMPROVE precompute the a^n_l up to a given limit and save them in a file to read?
"""
function fill_prefactors!(basis::structCB72Basis_type)

    lmax, nmax  = basis.lmax, basis.nmax
    G, rb       = basis.G, basis.rb
    tabPrefU, tabPrefD = basis.tabPrefU, basis.tabPrefD

    dimU = - sqrt(G / rb)                       # Potential basis element dimensional prefactor
    dimD = 1.0 / ( 2.0*pi*sqrt(G * (rb)^(3)) )  # Density basis element dimensional prefactor
    # Initialization
    for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        tabPrefU[1,n+1] = sqrt(2.0) * dimU
        tabPrefD[1,n+1] = sqrt(2.0) * dimD
    end 
    # Recurrence
    for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        for l=1:lmax # Loop over the harmonic indices. ATTENTION, harmonic index starts at l=0 (here 0 has been initialized)
            tabPrefU[l+1,n+1] = prefactors_recurrence_l_CB72(tabPrefU[l,n+1],l,n)
            tabPrefD[l+1,n+1] = prefactors_recurrence_l_CB72(tabPrefD[l,n+1],l,n)
        end
    end 

end


##################################################
# Computation of the potential basis elements
##################################################
"""
    rhoCB72(x)

returns the parameter -1 <= rho <= 1 for a given dimensionless radius x=r/rb.
"""
function rhoCB72(x::Float64)
    return (x*x - 1.0)/(x*x + 1.0) # Value of rho
end

"""
    U_init_l_CB72(x, l)

Gives the (n=0, l) adimensional potential basis element value at a given position x = r/rb,
xi^0_l(x) = prod(2k - 1, k=1 -> l) / [  (1 + x^2)^{l + 1/2}   ].

@IMPROVE other function to give the initialization for all l simultaneously through recurrence?
"""
function U_init_l_CB72(x::Float64,
                        l::Int64)
    return prod(1:2:(2*l-1)) / (sqrt(1.0 + x*x) * (1.0 + x*x)^(l)) 
end

"""
    U_recurrence_n_CB72(u^{n-2}, u^{n-1}, rho(x), l, n)

Gives the next adimensional potential basis element (recurrence over radial number n),
xi^n_l(x), given the 2 previous ones, xi^{n-1}_l(x) and xi^{n-2}_l(x).

Initializing the recurrence at xi^0_l given by the initialization from U_init_l_CB72,
one get the right next adimensional potential basis element from Clutton-Brock (1972).
"""
function U_recurrence_n_CB72(u0::Float64, u1::Float64, 
                                rho::Float64,
                                l::Int64, n::Int64)
    return ( 2.0 + ((2.0*l-1.0)/n) ) * rho * u1 - ( 1.0 + ((2.0*l-1.0)/n) ) * u0
end

"""
    tabUl!(basis::structCB72Basis_type, l, r[, forD])

For Clutton-Brock (1972) basis elements.
"""
function tabUl!(basis::structCB72Basis_type,
                    l::Int64,r::Float64,
                    forD::Bool=false)

    nmax        = basis.nmax
    rb          = basis.rb
    tabPrefU    = basis.tabPrefU
    tabl        = forD ? basis.tabDl : basis.tabUl

    x   = r/rb        # Dimensionless radius
    xl  = x^(l)
    rho = rhoCB72(x)  # Value of the rescaled parameter rho

    #####
    # Recurrence on adimensional subpart (xi) of the potential basis elements
    #####
    # Initialization
    u0, u1 = 0.0, U_init_l_CB72(x,l)  # n = -1, 0
    tabl[1] =  u1
    # Recurrence loop
    for n=1:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        u2 = U_recurrence_n_CB72(u0,u1,rho,l,n)
        tabl[n+1] = u2
        u0, u1 = u1, u2
    end

    #####
    # Multiplying by the appropriate prefactors (except if nopref=true, then return only xi^n_l(x))
    #####
    if !(forD)
        for n=0:nmax # Loop over the radial basis numbers. ATTENTION, index starts at n=0
            tabl[n+1] *= tabPrefU[l+1,n+1] * xl
        end
    end
end

"""
    getUln(basis::structCB72Basis_type, l, n, r[, forD])

For Clutton-Brock (1972) basis elements.
"""
function getUln(basis::structCB72Basis_type,
                    l::Int64,n::Int64,r::Float64,
                    forD::Bool=false)

    rb          = basis.rb
    tabPrefU    = basis.tabPrefU

    x   = r/rb        # Dimensionless radius
    xl  = x^(l)
    rho = rhoCB72(x)  # Value of the rescaled parameter rho

    #####
    # Recurrence on adimensional subpart (xi) of the potential basis elements
    #####
    # Initialization
    u0, u1 = 0.0, U_init_l_CB72(x,l)  # n = -1, 0
    # Recurrence loop
    for np=1:n # Loop over the radial basis numbers. ATTENTION, index starts at n=0
        u2 = U_recurrence_n_CB72(u0,u1,rho,l,np)
        u0, u1 = u1, u2
    end

    #####
    # Multiplying by the appropriate prefactor (except if nopref=true, then return only xi^n_l(x))
    #####
    if !(forD)
        u2 *= tabPrefU[l+1,n+1] * xl
    end

    return u2
end


##################################################
# Computation of the density basis elements
##################################################
"""
    tabDl!(basis::structCB72Basis_type, l, r)

For Clutton-Brock (1972) basis elements, using tabUl!.
"""
function tabDl!(basis::structCB72Basis_type,
                    l::Int64,r::Float64)
    #####
    # Compute potential basis elements at azimuthal number l+1 without prefactors 
    # Stored in tabD
    #####
    tabUl!(basis,l+1,r,true) 
    #####
    # Deduce density basis elements at azimuthal number l without prefactors 
    #####
    nmax        = basis.nmax
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
        tabDl[n+1] *= tabPrefD[l+1,n+1] * xl
    end
end

"""
    getDln(basis::structCB72Basis_type, l, n, r)

For Clutton-Brock (1972) basis elements, using getUln.
"""
function tabDl!(basis::structCB72Basis_type,
                    l::Int64,r::Float64)
    #####
    # Compute potential basis elements n-2 and n at azimuthal number l+1 without prefactors 
    #####
    u0, u2 = getUln(basis,l+1,n-2,r,true), getUln(basis,l+1,n,r,true)
    #####
    # # Deduce density basis elements at azimuthal number l without prefactors 
    #####
    Dln = u2 - u0 
    #####
    # Adding prefactors 
    #####
    rb          = basis.rb
    tabPrefD    = basis.tabPrefD

    x   = r/rb        # Dimensionless radius
    xl  = x^(l)
    Dln *= tabPrefD[l+1,n+1] * xl

    return Dln
end