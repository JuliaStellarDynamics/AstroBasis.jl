##################################################
# Construction of the radial basis elements, for empirical orthogonal functions
# from Weinberg (1999)
##################################################

"""
    EXPeofBasis_type

Radial basis elements derived by SLE solution, pre-tabulated
"""
struct EXPeofBasis_type <: AbstractAstroBasis

    name::String         # Basis name (default CB73)
    dimension::Int64     # Basis dimension (default 2)

    lmax::Int64         # Maximal harmonic/azimuthal index (starts at 0)
    nmax::Int64         # Maximal radial index (starts at 0)

    G::Float64       # Gravitational constant (MUST be 1.0)
    rb::Float64      # Radial extension (default 1.0)

    tabPrefU::Array{Float64,2}  # Potential prefactors array
    tabPrefD::Array{Float64,2}  # Density prefactors array

    tabUl::Array{Float64,1}     # Potential elements value array
    tabDl::Array{Float64,1}     # Density elements value array

end



#isocache = "/Volumes/External1/Isochrone/.slgrid_sph_isochrone"

# set the data path for the basis prefactor elements
data_path() = abspath(joinpath(@__DIR__, "tables", "slgrid_sph_isochrone"))
# rename for isochrone specifically

##################################################
function iso_potential(r::Float64)
    bc = 1.
    return -(1)/(bc+sqrt(bc*bc + r*r))
end

function iso_density(r::Float64)
    bc = 1.
    a = sqrt(r*r+bc*bc)
    return (3*(bc+a)*a*a - r*r*(bc+3*a))/(4*pi*((bc+a)^3)*a*a*a)
end
##################################################





##################################################
function read_cache(isocache::String)


    io = open(isocache,"r")

    seek(io,0)
    lmax = read(io,UInt32)
    nmax = read(io,UInt32)
    numr = read(io,UInt32)
    cmap = read(io,UInt32)

    rmin  = read(io,Float64)
    rmax  = read(io,Float64)
    scale = read(io,Float64)

    EVtable    = Array{Float64}(undef, (lmax+1,nmax))
    Rtable     = Array{Float64}(undef,(numr))
    XItable    = Array{Float64}(undef,(numr))
    EFtablePOT = Array{Float64}(undef, (lmax+1,nmax,numr))
    EFtableDEN = Array{Float64}(undef, (lmax+1,nmax,numr))



    #lmax2 = (lmax+1)*(lmax+1)
    #print(lmax2)

    xmin = (rmin/scale - 1.0)/(rmin/scale + 1.0)
    xmax = (rmax/scale - 1.0)/(rmax/scale + 1.0)
    dxi = (xmax-xmin)/(numr)


    for r=1:numr
        xi = xmin + r*dxi
        Rtable[r]  = ((1.0+xi)/(1.0-xi))*scale
        XItable[r] = xi
    end


    for l=1:lmax+1
        dummyl = read(io,UInt32)
        for n=1:nmax
            EVtable[l,n] = read(io,Float64)
        end
        for n=1:nmax
            for r=1:numr
                tmpval = read(io,Float64)
                EFtableDEN[l,n,r] = tmpval*sqrt.(EVtable[l,n])*  iso_density(Rtable[r])
                EFtablePOT[l,n,r] = tmpval/sqrt.(EVtable[l,n])*iso_potential(Rtable[r])
            end
        end
    end

    return Rtable,XItable,EFtableDEN,EFtablePOT
end
##################################################




##################################################
function read_cache_header(isocache::String)

    io = open(isocache,"r")

    seek(io,0)
    lmax = read(io,UInt32)
    nmax = read(io,UInt32)
    numr = read(io,UInt32)
    cmap = read(io,UInt32)

    rmin  = read(io,Float64)
    rmax  = read(io,Float64)
    scale = read(io,Float64)
    xmin = (rmin/scale - 1.0)/(rmin/scale + 1.0)
    xmax = (rmax/scale - 1.0)/(rmax/scale + 1.0)
    dxi = (xmax-xmin)/(numr)

    return lmax,nmax,scale,numr,xmin,xmax,dxi
end

##################################################


##################################################
function r_to_xi(r::Float64,scale::Float64)
    return (r/scale-1.0)/(r/scale+1.0)
end

function xi_to_r(xi::Float64,scale::Float64)
    return (1.0+xi)/(1.0 - xi) * scale;
end
##################################################

##################################################
function specify_index(x::Float64,XItable::Array{Float64,1})

    xmin = minimum(XItable)
    xmax = maximum(XItable)
    dxi  = XItable[2]-XItable[1]

    #print(xmin,xmax,dxi)
    indx = Int(floor((x-xmin)/dxi))

    if (indx<1)
        indx = 1
    end
    if (indx>numr-1)
        indx = numr-1;
    end

    x1 = (XItable[indx+1] - x)/(dxi);
    x2 = (x - XItable[indx])/(dxi);

    return indx,x1,x2
    # also consider x1,x2
end
##################################################


function fill_prefactors!(isocache::String)
    RTable,XItable,EFtableDEN,EFtablePOT = read_cache(isocache)
    lmax_eof,nmax_eof,scale,numr,xmin,xmax,dxi = read_cache_header(isocache)
end



##################################################
# Definition of the radial basis elements
# from SLE, Weinberg (1999)
# Arguments are:
# + l: harmonic index
# + np: radial index
# + r: radius
##################################################
function UlnpEOF(l::Int64,np::Int64,r::Float64)
    indx,x1,x2 = specify_index(r_to_xi(r,scale),XItable)
    return (x1*EFtablePOT[l+1,np,indx] + x2*EFtablePOT[l+1,np,indx+1])
end
#####
function DlnpEOF(l::Int64,np::Int64,r::Float64)
    indx,x1,x2 = specify_index(r_to_xi(r,scale),XItable)
    return (x1*EFtableDEN[l+1,np,indx] + x2*EFtableDEN[l+1,np,indx+1])
end

"""# Function that computes the values of Ulnp(r)
# for a given l and r,
# and for 1 <= np <= nradial
"""
function tabUlnpEOF!(l::Int64,r::Float64,
                       tabUlnp::Array{Float64,1})
    #####
    for np=1:nradial # Loop over the radial indices
        tabUlnp[np] = UlnpEOF(l,np,r) # Filling the value for np
    end
end
