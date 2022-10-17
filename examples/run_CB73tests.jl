"""
AstroBasis Clutton-Brock (1973) test file


"""

import AstroBasis

# for robust benchmarking; could skip this
using BenchmarkTools

# set up the basis
rb,G = 1.,1.

CB73 = AstroBasis.CB73BasisCreate(lmax=6, nmax=20,G=1., rb=1.)
AstroBasis.getUln(CB73,1,1,0.1)

AstroBasis.tabUl!(CB73,0,0.1)
println(CB73.tabUl)

AstroBasis.tabDl!(CB73,0,0.1)

# run some timing tests
function EvalTime(l::Int64,nr::Int64,basis::AstroBasis.BasisType)
    for x=1:nr
        r = x*0.001
        AstroBasis.tabUl!(basis,0,0.1)
    end
end

@btime EvalTime(0,1000,CB73)

# using a different basis is very straightforward, e.g.
#HO92 = AstroBasis.HernquistBasisCreate(lmax=6, nmax=20,G=1., rb=1.)
#@btime EvalTime(0,1000,HO92)
