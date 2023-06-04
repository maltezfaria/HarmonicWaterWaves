using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 5,10
h = 0.4
q = 3
θ = π/10
d = 2

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

# create piercing obstacle
r = 1.0 # radius of obstacle
δ = 3*r
x₀  = SVector(-δ/2,0.0) # center of obstacle
obs = let r = r, x₀ = x₀
    WPB.ParametricEntity(0,π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₀
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

x₁  = SVector(δ/2,0.0) # center of obstacle
obs = let r = r, x₁ = x₁
    WPB.ParametricEntity(0,π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₁
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

# free surface
HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,-r-δ/2)
HarmonicWaterWaves.add_freesurface!(tank,-δ/2 + r, δ/2 - r)
HarmonicWaterWaves.add_freesurface!(tank,δ/2+r,a)
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

idxs_free   = WPB.dom2qtags(tank.quad,WW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(tank.quad,WW.bottom(tank))
idxs_obs   = WPB.dom2qtags(tank.quad,WW.obstacles(tank))

# solve
HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
F = HarmonicWaterWaves.solve_eigenvalues(tank)
λ = filter(x->!isnan(x),F.values)
scatter(λ,label="")
