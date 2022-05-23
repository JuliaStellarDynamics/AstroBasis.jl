"""
Radial basis creation requisites

- Structure with the following attributes :
    - name
    - dimension
    - lmax
    - nmax
    - rb
    - G

    - tabPrefU
    - tabPrefD
    - tabUl
    - tabDl

@WARNING: Arrays index are offsetted as basis usually start at l=0, n=0.

- Adapted following functions :
    - fill_prefactors!
    - getUln
    - getDln
    - tabUl!
    - tabDl!

"""


"""fill_prefactors!(basis)

Fill the density and potential prefactors arrays from 'basis'

@WARNING: this may only be run once per session, as it will create constants.
"""
function fill_prefactors!(basis)
    # ... [implementation sold separately] ...
end


"""getUln(basis,l,n,r)

Definition of the potential radial basis elements (l,n) at position r.

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function getUln(basis,l,n,r)
    # ... [implementation sold separately] ...
end


"""getDln(basis,l,n,r)

Definition of the density radial basis elements (l,n) at position r.

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function getUln(basis,l,n,r)
    # ... [implementation sold separately] ...
end


"""tabUl!(basis,l)

Compute potential table for a given l and r, and for 0 <= n <= nmax

To avoid repeated memory allocations, this overwrites the 'basis.tabUl' table

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function tabUl!(basis,l,r)
    # ... [implementation sold separately] ...
end

"""tabDl!(basis,l)

Compute density table for a given l and r, and for 0 <= n <= nmax

To avoid repeated memory allocations, this overwrites the 'basis.tabDl' table

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function tabDl!(basis,l,r)
    # ... [implementation sold separately] ...
end