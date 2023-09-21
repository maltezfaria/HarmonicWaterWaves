# HarmonicWaterWaves

**Boundary integral equation methods for solving time-harmonic water waves in
Julia**

## Installing

To add this package, type the following code on a Julia REPL:

```julia
using Pkg
Pkg.add("https://github.com/maltezfaria/HarmonicWaterWaves.git")
```

## Usage

Currently, this package is only intended as a reproducibility repository for the
paper [*Complex-scaled boundary integral equation for time-harmonic water
waves*](paper/tex/water-waves-pml.pdf). Run the following code on a terminal
from the root of this repository to regenerate all the figures presented in the
paper:

```bash
julia --project=paper paper/makefigures.jl
```

This will repopulate the [figures](figure) and [animations](animations) folders.
If you have a working `LaTex` installation, you can recompile the
[`tex/water-waves-pml.tex`](tex/water-waves-pml.tex) file to produce the final
`.pdf`.
