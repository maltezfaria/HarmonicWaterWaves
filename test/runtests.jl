using HarmonicWaterWaves
using Test
using StaticArrays
using HCubature
using LinearAlgebra
using Nystrom
using Nystrom: ParametricSurfaces
using Plots
import WavePropBase as WPB

plotlyjs()

include(joinpath(HarmonicWaterWaves.PROJECT_ROOT,"test/exact_solution.jl"))

a, l = 10, 20
h = 0.05
q = 3

p    = HarmonicWaterWaves.Parameters(frequency=√π,gravity=1,θ=π/4)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.add_free_surface!(tank,a)
HarmonicWaterWaves.add_pml_layer!(tank,l)
obs = boundary(ParametricSurfaces.Disk(;radius=0.5,center=(0.0,-1.0),normal=:inside))[1]
HarmonicWaterWaves.add_obstacle!(tank,obs)

plot(tank)

HarmonicWaterWaves.discretize!(tank;meshsize=h,order=q)

HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)

function f(dof)
    x = dof.coords
    if x[2] == 0
        return 0.0
    else
        return 1.0
    end
end

ϕ = HarmonicWaterWaves.solve(tank,f)

Γ = HarmonicWaterWaves.freesurface(tank)

quad   = tank.quad
xx   = [dof.coords[1] for dof in quad.dofs]

fig = plot()
for d in Γ
    idxs = Nystrom.dom2dof(quad,Domain(d))
    plot!(fig,xx[idxs],real(ϕ.vals[idxs]),label="",color=:red)
end
fig
