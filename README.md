# HarmonicWaterWaves

*Boundary integral equation methods for solving time-harmonic water waves in
Julia*

## Usage

Currently, this package is only intended as a reproducibility repository for the
paper [*Complex-scaled boundary integral equation for time-harmonic water
waves*](paper/tex/water-waves-pml.pdf). It requires installing [`julia`](https://julialang.org/downloads/),
version `1.9`.

Run the following code on a terminal to download the repository and regenerate
all the figures presented in the paper (you must have `git` installed):

```bash
git clone https://github.com/maltezfaria/HarmonicWaterWaves.git
cd ./HarmonicWaterWaves
julia --project=paper -e 'using Pkg; Pkg.instantiate(); include("paper/makefigures.jl")'
```

**:warning: Because of some (heavy) plotting dependencies that need to be downloaded and
precompiled, the lines above may take some time the first time you run it.**

This will repopulate the [figures](paper/figure) and [animations](paper/animations) folders.
If you have a working `LaTex` installation, you can recompile the
[`paper/tex/water-waves-pml.tex`](paper/tex/water-waves-pml.tex) file to produce the final
`.pdf`. Below is an example of the results produced (see the paper for more
details).

|      Obstacle scattering       |       Step scattering       |
| :----------------------------: | :-------------------------: |
| ![jelly](jellyfish_fields.gif) | ![step](step_animation.gif) |
