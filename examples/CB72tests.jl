
import AstroBasis
using BenchmarkTools

rb,G = 1.,1.
basis = AstroBasis.CB72BasisCreate(lmax=2, nmax=100,G=1., rb=1.)

r = 1.5
l = 2

# @benchmark AstroBasis.tabUl!(basis,l,r)