"""
AstroBasis Clutton-Brock (1973) test file

@ATTENTION: to check pre-n-shift values, try
git checkout 2cbfb7d8d1bbe9c570331a282bf161297e7d37cd


"""

import AstroBasis

# set up the basis
rb,G = 1.,1.

CB73 = AstroBasis.CB73Basis_create(lmax=6, nmax=20,G=1., rb=1.)
AstroBasis.fill_prefactors!(CB73)
AstroBasis.getUln(CB73,1,1,0.1)


AstroBasis.tabUl!(CB73,0,0.1)
println(CB73.tabUl)

AstroBasis.tabDl!(CB73,0,0.1)

# run some timing tests
function eval_time(l::Int64,nr::Int64,basis::AstroBasis.structCB73Basis_type)
    for x=1:nr
        r = x*0.001
        AstroBasis.tabUl!(basis,0,0.1)
    end
end

@time eval_time(0,1000,CB73)
