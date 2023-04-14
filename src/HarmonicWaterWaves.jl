module HarmonicWaterWaves

const PROJECT_ROOT =  pkgdir(HarmonicWaterWaves)

using LinearAlgebra
using Roots
using RecipesBase
using StaticArrays

# TODO: make a PR to SpecialFunctions to support complex expinti?
import SpecialFunctions
SpecialFunctions.expinti(z::Complex) = -SpecialFunctions.expint(-z) - sign(angle(-z))*im*π

import WavePropBase:
    Domain,
    coords,
    svector,
    jacobian,
    AbstractPDE,
    SingleLayerKernel,
    DoubleLayerKernel,
    default_kernel_eltype,
    default_density_eltype,
    ambient_dimension,
    notimplemented,
    normal,
    line,
    Point2D,
    meshgen,
    NystromMesh,
    Laplace,
    SingleLayerOperator,
    DoubleLayerOperator,
    hcubature_correction,
    dom2qtags,
    NystromDensity



include("pml.jl")
include("parameters.jl")
include("wavetank.jl")

# export

end
