
# AstroBasis.jl

[![image](https://github.com/JuliaStellarDynamics/AstroBasis.jl/actions/workflows/documentation.yml/badge.svg?branch=documentation)](https://juliastellardynamics.github.io/AstroBasis.jl/)

`AstroBasis.jl` is a package written in Julia to compute basis function evaluations, primarily for astronomical models.

---
## Installation

Install Julia by following the instructions at [julialang.org/downloads/](https://julialang.org/downloads/).

To invoke Julia in the Terminal, you need to make sure that the `julia` command-line program is in your `PATH`. 
See [here](https://julialang.org/downloads/platform/#optional_add_julia_to_path) for detailed instructions.

*If you do not want to install Julia but want to test the library, you can use this [Google colab notebook](https://colab.research.google.com/drive/1g5AD8zzwyqmufqVdYEzkdi5hdifu-z2S?usp=sharing).
However, Google colab is not primarly made to run Julia code. 
It will then need to be installed on the remote machine which can take a few minutes.
This notebook is not maintained as a priority. We would recommand you install Julia on your machine to test the library locally.*

Once Julia installed, you need to install the `AstroBasis.jl` library. You can either 
* add it to your global Julia environment by simply running:
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
*Note that if you use this second install option you will always need to run codes in the project context by adding the option `--project=/path/to/AstroBasis.jl` after `julia`. 
The library will not be accessible in your global Julia context.*

Then, in your program, if you want use `AstroBasis` function, write `import AstroBasis`.

---
## Quick use test

An introduction example is given in `example/test_AstroBasis.jl`.

If you installed the library using the first (global) install option, just download this example [file](https://github.com/JuliaStellarDynamics/AstroBasis.jl/blob/main/examples/test_AstroBasis.jl) from the github repository.

Run the code with the following command:
```
$ julia /path/to/test_AstroBasis.jl
```
*Do not forget the option `--project=/path/to/AstroBasis.jl` after `julia` if you installed it the local way.*

This example will first install some required libraries (`Plots`, `LaTeXStrings`). These installations might take up to a minute when first called.
If you want to modify the parameters in the example file and run it again, you can comment the corresponding lines to avoid installation checks.

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
