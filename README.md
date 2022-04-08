
# AstroBasis.jl

`AstroBasis.jl` is a package written in Julia to compute basis function evaluations, primarily for astronomical models.

-----------------------------

## Quick activate

In the main directory where you the package lives, enter the Julia environment (`julia`), then the package manager (`]`), then activate (`activate .`). To be extra safe, you can `resolve` to check for updates. Then return to the Julia interpreter (`[backspace]`): you are good to go with the latest version of the package! Import the exports by typing `using AstroBasis` into the Julia interpreter. You may also need to download some packages if you are using a new Julia interpreter: try `using(Pkg);Pkg.instantiate()`. If you want to access specific elements listed below, I recommend `import AstroBasis`.

As `AstroBasis` is (currently) unregistered, if you would like to add it to your Julia registry, read [here](https://pkgdocs.julialang.org/v1/managing-packages/#Adding-unregistered-packages). Short version: when in the package manager, `add "git@github.com:michael-petersen/JuliaAstroBasis.git"`. If you are getting an error about git keys, you will need to register your private key using the julia shell prompt (access with `;`), and then pointing at your private key: `ssh-add ~/.ssh/id_rsa`.

-----------------------------

## Clutton-Brock 3d basis functions

`UlnpCB73(l,n,r,tabUln)` will calculate the Clutton-Brock basis at location r for harmonic order l and radial order n. Requires tabulated potential prefactors, tabUln.

`read_and_fill_prefactors(lmax,nmax)` will load prefactors for the Clutton-Brock basis.

-----------------------------

## Author

Mike Petersen -  @michael-petersen - petersen@iap.fr
