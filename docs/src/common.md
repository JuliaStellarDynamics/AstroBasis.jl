# Functions common to all bases

To seamlessly allow for switching between different bases, the bases are all organised around common functions. 

## Common functions

The common functions are detailed in `src/Basisdoc.jl`. Any instantiated basis will result in a structure, with attributes `name`, `dimension`, `lmax`, `nmax`, `rb` (the basis scale), and `G` (the gravitational constant; unity by default). The following functions are also defined:

`getUln(basis,l,n,r)` will calculate the basis potential at location r for harmonic order l and radial order n.

`getDln(basis,l,n,r)` will calculate the basis density at location r for harmonic order l and radial order n.

Overloaded function `tabUl!(basis,l)` will fill the structure attribute `tabUl` table with radial potential values from the basis up to `nmax`.

Overloaded function `tabDl!(basis,l)` will fill the structure attribute `tabUl` table with radial density values from the basis up to `nmax`.


## Accessing potential functions
```@docs
AstroBasis.getUln
```

```@docs
AstroBasis.tabUl!
```

## Accessing density functions

```@docs
AstroBasis.getDln
```

```@docs
AstroBasis.tabDl!
```
