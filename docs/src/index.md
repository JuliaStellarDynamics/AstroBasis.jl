# AstroBasis.jl

*Bi-orthogonal bases for galactic dynamics in Julia*

---
## Installation

For installation and first use test, follow the instructions on the github [README](https://github.com/JuliaStellarDynamics/AstroBasis.jl).

---
## Basis functions

`AstroBasis` currently offers four different types of bases.
### 3d
1. Clutton-Brock (1973), which is a match to the Plummer (1911) profile.
2. Hernquist & Ostriker (1992), which is a match to the Hernquist (1990) profile.

### 2d
1. Clutton-Brock (1972), which is a 2d match to a Plummer profile.
2. Kalnajs (1976), which is a spiral basis.

---

## Common functions

The common functions are detailed in `src/Basisdoc.jl`. Any instantiated basis will result in a structure, with attributes `name`, `dimension`, `lmax`, `nmax`, `rb` (the basis scale), and `G` (the gravitational constant; unity by default). The following functions are also defined:

`getUln(basis,l,n,r)` will calculate the basis potential at location r for harmonic order l and radial order n.

`getDln(basis,l,n,r)` will calculate the basis density at location r for harmonic order l and radial order n.

Overloaded function `tabUl!(basis,l)` will fill the structure attribute `tabUl` table with radial potential values from the basis up to `nmax`.

Overloaded function `tabDl!(basis,l)` will fill the structure attribute `tabUl` table with radial density values from the basis up to `nmax`.

`read_and_fill_prefactors(lmax,nmax)` will load prefactors for the basis. This is called when the basis is instantiated, and is unlikely to be needed to be called externally.
