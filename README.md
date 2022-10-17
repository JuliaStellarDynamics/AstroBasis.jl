
# AstroBasis.jl

`AstroBasis.jl` is a package written in Julia to compute basis function evaluations, primarily for astronomical models.

-----------------------------

### Quick activate

`AstroBasis` is (currently) unregistered, and as such if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Shortcut version: after cloning the repository, navigate to the top directory. Start julia (`julia`), then the package manager (`]`), then register the package in development mode (`dev .`).

Then, in your program, if you want to access specific elements listed below, use `import AstroBasis`.

-----------------------------

### Basis functions

`AstroBasis` currently offers four different types of bases.
#### 3d
1. Clutton-Brock (1973), which is a match to the Plummer (1911) profile.
2. Hernquist & Ostriker (1992), which is a match to the Hernquist (1990) profile.

#### 2d
1. Clutton-Brock (1972), which is a 2d match to a Plummer profile.
2. Kalnajs (1976), which is a spiral basis.

-----------------------------

### Common functions

The common functions are detailed in `src/Basisdoc.jl`. Any instantiated basis will result in a structure, with attributes `name`, `dimension`, `lmax`, `nmax`, `rb` (the basis scale), and `G` (the gravitational constant; unity by default). The following functions are also defined:

`getUln(basis,l,n,r)` will calculate the basis potential at location r for harmonic order l and radial order n.

`getDln(basis,l,n,r)` will calculate the basis density at location r for harmonic order l and radial order n.

Overloaded function `tabUl!(basis,l)` will fill the structure attribute `tabUl` table with radial potential values from the basis up to `nmax`.

Overloaded function `tabDl!(basis,l)` will fill the structure attribute `tabUl` table with radial density values from the basis up to `nmax`.

`read_and_fill_prefactors(lmax,nmax)` will load prefactors for the basis. This is called when the basis is instantiated, and is unlikely to be needed to be called externally.

-----------------------------

### Examples

See `examples/run_CB73tests.jl` for an example script using the Clutton-Brock 3d basis and a timing test. The example is generic to all bases: one can simply swap the initialiser.

-----------------------------

### Authors

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk, petersen@iap.fr

Mathieu Roule -  @MathieuRoule - roule@iap.fr
