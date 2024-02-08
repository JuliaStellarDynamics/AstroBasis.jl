
# AstroBasis.jl

[![image](https://github.com/JuliaStellarDynamics/AstroBasis.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/AstroBasis.jl/)

`AstroBasis` is a package written in Julia to compute basis function evaluations, primarily for astronomical models.

---
## Quick installation and use test

Install Julia by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/).

To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

Once Julia installed, clone the `AstroBasis.jl` library and precompile it by running:
```
git clone https://github.com/JuliaStellarDynamics/AstroBasis.jl.git
cd AstroBasis.jl
julia --project=. -e 'using Pkg; Pkg.precompile()'
```

An introduction example is given in `example/test_AstroBasis.jl`.
Run the code with the following command:
```
julia --project=. examples/test_AstroBasis.jl
```
*Make sure you still are at the* `/path/to/AstroBasis.jl` *location.*

This example will first install some required libraries (`Plots`, `LaTeXStrings`). These installations might take up to a minute when first called.

The resulting plot will be created in the same folder as the test code under the name `CluttonBrock73.png`.

![`Clutton-Brock (1973)`](examples/CluttonBrock73_original.png)


### Interactive notebook

If you prefer interactive Jupyter notebooks, you will need to install `IJulia` following these [instructions](https://github.com/JuliaLang/IJulia.jl).

The interactive introduction example is then given in `example/test_AstroBasis.ipynb`.

---
### Without installing Julia

*If you do not want to install Julia but want to test the library, you can use this [Google colab notebook](https://colab.research.google.com/drive/1g5AD8zzwyqmufqVdYEzkdi5hdifu-z2S?usp=sharing).
However, Google colab is not primarly made to run Julia code. 
It will then need to be installed on the remote machine which can take a few minutes.
This notebook is not maintained as a priority. We would recommand you install Julia on your machine to test the library locally.*

---
## More recurrent usage

If you followed the quick installation and use test, you have installed `ÀstroBasis` and ran the example in a local environment (or project) associated to the clone's folder.
The library then won't be accessible outside this local scope (hence the need for the `--project=.` option).
If you want to use the library functions more regularly, we recommand you bring it to the global scope by either
* linking the clone to the global scope by running
```
julia -e 'using Pkg; Pkg.add(url="/localpath/to/AstroBasis.jl")'
```
*This link might be affected by any change in the clone's location.*

or 
* installing the library directly in Julia's global scope (it will clone the repository somewhere in Julia's folder `.julia/packages`) by running
```
julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/AstroBasis.jl.git")'
```
*You can then delete your "local" clone.*

You will then be able to use `AstroBasis` functions in any of your scripts by writing `import AstroBasis`.

---
## Documentation and usage

To get more familiar with the content of the library and start and design your own use case, you may want to visit the [documentation](https://juliastellardynamics.github.io/AstroBasis.jl/).


---
## Authors

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk, petersen@iap.fr

Mathieu Roule -  @MathieuRoule - roule@iap.fr
