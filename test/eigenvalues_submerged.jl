using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 0.55,10
h = 0.2
q = 4
θ = π/4
d = 2

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

# create piercing obstacle
r = 0.5 # radius of obstacle
x₀  = SVector(0,-1.0) # center of obstacle
obs = let r = r, x₀ = x₀
    WPB.ParametricEntity(0,2π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₀
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

# free surface
HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

# bottom
HarmonicWaterWaves.set_depth!(tank,d)
HarmonicWaterWaves.add_bottom!(tank,-a-l,-a)
HarmonicWaterWaves.add_bottom!(tank,-a,a)
HarmonicWaterWaves.add_bottom!(tank,a,a+l)

# add pml
HarmonicWaterWaves.add_orthogonal_pml!(tank;a,θ)

# create mesh
HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

# solve
HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
F = HarmonicWaterWaves.solve_eigenvalues(tank)

##
scatter(F.values,label="")
