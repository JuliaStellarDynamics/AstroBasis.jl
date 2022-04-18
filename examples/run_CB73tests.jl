"""
AstroBasis Clutton-Brock (1973) test file

@IMPROVE: test cumsum on arrays versus loop for timing
"""

import AstroBasis

# set up the basis
rb,G = 1.,1.
UF,DF = AstroBasis.read_and_fill_prefactors(6,20,rb,G)

# return a single basis function value
println(AstroBasis.UlnpCB73(1,2,0.1,UF))

# initialise an array for holding the basis values: memory efficient!
tabUlnp = zeros(20)
AstroBasis.tabUlnpCB73!(0,0.1,tabUlnp,20,UF,rb)

# there could also be a version that returns the array, just as a backup for one-offs.

# run some timing tests
function eval_time(l::Int64,nr::Int64,
                   tabUlnp::Array{Float64,1},
                   nradial::Int64,
                   tabPrefCB73_Ulnp::Matrix{Float64},
                   rb::Float64=1.)
    for x=1:nr
        r = x*0.001
        AstroBasis.tabUlnpCB73!(0,r,tabUlnp,20,tabPrefCB73_Ulnp)
    end
end

@time eval_time(0,1000,tabUlnp,20,UF)
