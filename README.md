
# AstroBasis.jl

[![image](https://github.com/JuliaStellarDynamics/AstroBasis.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/AstroBasis.jl/)

`AstroBasis` is a package written in Julia to compute basis function evaluations, primarily for astronomical models.

---
## Quick installation and use test

Install Julia by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/).

To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

We will now proceed to the installation of the `AstroBasis` library.

**Caution:** If you wish to exercise extra care to prevent any interference with other installed libraries (e.g., avoiding updates), please refer to the "Notes on working with environments" section below.

Install the `AstroBasis` library by running
```
julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/AstroBasis.jl.git")'
```

An introduction example is given in `examples/test_AstroBasis.jl`. Download the file by running:
```
wget https://raw.githubusercontent.com/JuliaStellarDynamics/AstroBasis.jl/main/examples/test_AstroBasis.jl
```
Then run the code with the following commands
```
julia test_AstroBasis.jl
```

This example will first install some required libraries (`Plots`, `LaTeXStrings`) and their dependencies. These installations might take up to 4 minutes.

The resulting plot will be created in the same folder as the test code under the name `CluttonBrock73.png`.

![`Clutton-Brock (1973)`](examples/CluttonBrock73_original.png)

---
### Note on working with environments

By default, packages are added to the default environment at ~/.julia/environments/v1.#.
It is however easy to create other, independent, projects.
If you want to install the `AstroBasis` package in a different/test environment, first create a folder to host the environment files (Project.toml and Manifest.toml which will be created later on).
Then, for every command line invoking Julia, use `julia --project=/path/to/my_env` instead of `julia` alone.

*Note that packages will always be cloned in ~/.julia/packages but only accessible in your project's context.* 
*A procedure to fully uninstall the package is described at the end of this readme.*

---
### Interactive notebook

If you prefer interactive Jupyter notebooks, you will need to install `IJulia` following these [instructions](https://github.com/JuliaLang/IJulia.jl).

An interactive introduction example is given in [`examples/test_AstroBasis.ipynb`](examples/test_AstroBasis.ipynb).

---
### Without installing Julia

*If you do not want to install Julia but want to test the library, you can use this [Google colab notebook](https://colab.research.google.com/drive/1g5AD8zzwyqmufqVdYEzkdi5hdifu-z2S?usp=sharing).
However, Google colab is not primarly made to run Julia code. 
It will then need to be installed on the remote machine which can take a few minutes.
This notebook is not maintained as a priority. We would recommand you install Julia on your machine to test the library locally.*

---
## Documentation and usage

To get more familiar with the content of the library and start and design your own use case, you may want to visit the [documentation](https://juliastellardynamics.github.io/AstroBasis.jl/).

---
## Uninstall

First start by removing the package from the environment by running
```
julia -e 'using Pkg; Pkg.rm("AstroBasis");'
```

Following the same syntax, you can also remove the `Plots` and `LaTeXString` packages installed for the example if you want to. 

If you worked in a test environment (that you do not want to keep) you can also simply erase the folder using `rm -r /path/to/my_env`.

Then to fully erase the package (installed in ~/.julia), run
```
julia -e 'using Pkg; using Dates; Pkg.gc(collect_delay=Day(0));'
```
<sup><sub>*No need for the* `--project=/path/to/my_env` *option here !*</sub></sup>

It will erase all the packages which are not known in any of your "active" (i.e., for which the Manifest.toml file is reachable) project/environments, in particular `AstroBasis`.

---
## Authors

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk, petersen@iap.fr

Mathieu Roule -  @MathieuRoule - roule@iap.fr
