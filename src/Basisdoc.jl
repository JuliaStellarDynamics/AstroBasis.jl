"""
    Radial basis creation requisites

- Substructure of the abstract type AbstractAstroBasis
with the following attributes :
    - name
    - dimension
    - lmax
    - nradial
    - rb
    - G
    - other usefull parameters

    - tabUl
    - tabDl
    - other usefull arrays 

@WARNING: mandatory fields → name, dimension, lmax, nradial, rb, G, tabUl, tabDl

Abstract type :
    struct MyBasisStruct <: AbstractAstroBasis
         
        myparameters...

        myarrays...
    end

- Adapted following methods :
    - structure constructor
    - getUln
    - getDln
    - tabUl!
    - tabDl!

"""


"""
    getUln(basis,l,n,r)

Definition of the potential radial basis elements (l,n) at position r.

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function getUln(basis,l,n,r)
    # ... [implementation sold separately] ...
end


"""
    getDln(basis,l,n,r)

Definition of the density radial basis elements (l,n) at position r.

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function getDln(basis,l,n,r)
    # ... [implementation sold separately] ...
end


"""
    tabUl!(basis,l)

Compute potential table for a given l and r, and for 0 <= n <= nmax

To avoid repeated memory allocations, this overwrites the 'basis.tabUl' table

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function tabUl!(basis,l,r)
    # ... [implementation sold separately] ...
end

"""
    tabDl!(basis,l)

Compute density table for a given l and r, and for 0 <= n <= nmax

To avoid repeated memory allocations, this overwrites the 'basis.tabDl' table

@WARNING:: l=0 is index 1, n=0 is index 1.
"""
function tabDl!(basis,l,r)
    # ... [implementation sold separately] ...
end
