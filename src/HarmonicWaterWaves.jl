module HarmonicWaterWaves

const PROJECT_ROOT =  pkgdir(HarmonicWaterWaves)

using LinearAlgebra
using Roots
using RecipesBase
using StaticArrays
using LinearMaps
using SparseArrays
using IterativeSolvers

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
    Point3D,
    ParametricEntity,
    meshgen,
    NystromMesh,
    Laplace,
    dom2qtags,
    SingleLayerOperator,
    DoubleLayerOperator,
    hcubature_correction,
    NystromDensity,
    IntegralPotential

include("pml.jl")
include("parameters.jl")
include("wavetank.jl")

# export

end
