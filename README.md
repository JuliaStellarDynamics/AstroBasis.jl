
# AstroBasis.jl

[![image](https://github.com/JuliaStellarDynamics/AstroBasis.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/AstroBasis.jl/)

`AstroBasis.jl` is a package written in Julia to compute basis function evaluations, primarily for astronomical models.

---
## Installation

If you are new to `julia`, install the latest version by running this in your terminal:
```
$ curl -fsSL https://install.julialang.org | sh
```
<font size="1"> If you are on Windows or struggle with installation, please visit [this website](https://julialang.org/downloads/). </font>

Once `julia` installed, you need to install this library. You can either 
* add it to your global julia environment by simply running:
```
$ julia -e 'using Pkg; Pkg.add(url="https://github.com/JuliaStellarDynamics/OrbitalElements.jl.git")'
```
or 
* clone the repository wherever you want and create a local environment (or project) by running:
```
$ git clone https://github.com/JuliaStellarDynamics/AstroBasis.jl.git
$ cd AstroBasis.jl
$ julia --project=. -e 'using Pkg; Pkg.precompile()'
```
<font size="1"> Note that if you use this second install option you will always need to run codes in the project context by adding the option `--project=/path/to/AstroBasis.jl` after `julia`. The library will not be accessible in your global julia context.</font>

Then, in your program, if you want use `AstroBasis` function, write `import AstroBasis`.

---
## Quick use test

An introduction example is given in `example/test_AstroBasis.jl`.

If you installed the library using the first (global) install option, just download this example [file](https://github.com/JuliaStellarDynamics/AstroBasis.jl/blob/main/examples/CB72tests.jl) from the github repository.

Run the code with the following command:
```
$ julia /path/to/test_AstroBasis.jl
```
<font size="1"> Do not forget the option `--project=/path/to/AstroBasis.jl` after `julia` if you installed it the local way.</font>

This example will first install some required libraries (`Plots`, `LaTeXStrings`). These installations might take a few minutes when first called.

The resulting plot will be created in the same folder as the test code under the name `CluttonBrock73.png`.

![`Clutton-Brock (1973)`](examples/CluttonBrock73.png)

### Interactive notebook

If you prefer interactive Jupyter notebooks, you will need to install `IJulia` following these [instructions](https://github.com/JuliaLang/IJulia.jl).

The interactive introduction example is then given in `example/test_AstroBasis.ipynb`.

---
## Documentation and usage

To get more familiar with the content of the library and start and design your own use case, you may want to visit the [documentation](https://juliastellardynamics.github.io/AstroBasis.jl/).


---
## Authors

Mike Petersen -  @michael-petersen - michael.petersen@roe.ac.uk, petersen@iap.fr

Mathieu Roule -  @MathieuRoule - roule@iap.fr
