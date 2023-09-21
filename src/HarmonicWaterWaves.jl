module HarmonicWaterWaves

const PROJECT_ROOT =  pkgdir(HarmonicWaterWaves)

using LinearAlgebra
using Roots
using StaticArrays
using LinearMaps
using SparseArrays
using IterativeSolvers
using CairoMakie
using Unitful

import WavePropBase as WPB

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
    hcubature_correction,
    SingleLayerOperator,
    DoubleLayerOperator,
    NystromDensity,
    IntegralPotential

include("pml.jl")
include("wavetank.jl")
include("geometries.jl")
include("makietheme.jl")

# export

end
