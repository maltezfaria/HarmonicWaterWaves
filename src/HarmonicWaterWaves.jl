module HarmonicWaterWaves

const PROJECT_ROOT =  pkgdir(HarmonicWaterWaves)

using LinearAlgebra
using Roots
using RecipesBase
using Nystrom
using StaticArrays
import Nystrom: ParametricSurfaces

# TODO: make a PR to SpecialFunctions
import SpecialFunctions

SpecialFunctions.expinti(z::Complex) = -SpecialFunctions.expint(-z) - sign(angle(-z))*im*π

import WavePropBase:
    Domain,
    Point1D,
    Point2D,
    AbstractEntity,
    HyperRectangle,
    jacobian,
    domain,
    ambient_dimension,
    normal,
    coords,
    mesh

include("pml.jl")
include("parameters.jl")
include("wavetank.jl")

# export

end
