# HarmonicWaterWaves

*Boundary integral equation methods for solving time-harmonic water waves in
Julia*

## Usage

Currently, this package is only intended as a reproducibility repository for the
paper [*Complex-scaled boundary integral equation for time-harmonic water
waves*](https://arxiv.org/pdf/2310.04127.pdf). It requires installing [`julia`](https://julialang.org/downloads/),
version `1.9`.

Run the following code on a terminal to download the repository and regenerate
all the figures presented in the paper (you must have `git` installed):

```bash
git clone https://github.com/maltezfaria/HarmonicWaterWaves.git
cd HarmonicWaterWaves
git checkout v0.1
julia paper/makefigures.jl
```

**:warning: Because of some (heavy) plotting dependencies that need to be downloaded and
precompiled, the lines above may take some time the first time you run it.**

This will populate the [figures](paper/figure) and [animations](paper/animations) folders.
Below is an example of the results produced (see the paper for more
details).

|      Obstacle scattering       |       Step scattering       |
| :----------------------------: | :-------------------------: |
| ![jelly](jellyfish_fields.gif) | ![step](step_animation.gif) |
