using HarmonicWaterWaves
using Test
using StaticArrays
using HCubature
using LinearAlgebra
using Nystrom
using Nystrom: ParametricSurfaces
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 2, 5
h = 0.2
q = 3

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1,θ=π/4)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.add_free_surface!(tank,a)
HarmonicWaterWaves.add_pml_layer!(tank,l)
x₀  = SVector(1,-1.0) # center of obstacle
# xₛ  = SVector(0,-1.0) # center of source
xₛ  = x₀
obs = boundary(ParametricSurfaces.Disk(;radius=0.5,center=x₀,normal=:inside))[1]
HarmonicWaterWaves.add_obstacle!(tank,obs)

HarmonicWaterWaves.discretize!(tank;meshsize=h,order=q)

k = WW.impedance(tank)
op  = WW.InfiniteDepthWaterWaves(;k)
G   = SingleLayerKernel(op)
dG  = DoubleLayerKernel(op)
ϕe   = source -> G(xₛ,source)
dϕe  = source -> dG(xₛ,source)

## Solve with proper Green function
Γs         = WW.obstacles(tank)
quad_obs   = NystromMesh(view(tank.mesh,Γs);order=q)
op         = WW.InfiniteDepthWaterWaves(;k)
S          = SingleLayerOperator(op,quad_obs) |> Nystrom.assemble_gk
D          = DoubleLayerOperator(op,quad_obs) |> Nystrom.assemble_gk

f          = [dϕe(dof) for dof in quad_obs.dofs]

ϕ_green = (I/2 + D)\(S*f)

xx   = [atan(dof.coords[2]-x₀[2],dof.coords[1]-x₀[1]) for dof in quad_obs.dofs]
yy   = [ϕe(dof) for dof in quad_obs.dofs]

@info norm(yy-ϕ_green)

scatter(xx,real(yy),label="exact",m=:cross,ms=10)
scatter!(xx,real(ϕ_green),label="green",m=:circle)
##

## Solve with PML Green function
# plot(tank)

quad = tank.quad
dofs = quad.dofs

HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)

f = [dof.coords[2] == 0 ? dϕe(dof) - k*ϕe(dof) : dϕe(dof) for dof in dofs]

ϕ_pml = HarmonicWaterWaves.solve(tank,f)

idxs = Nystrom.dom2dof(quad,Γs)

@info norm(yy-ϕ_pml[idxs])

scatter!(xx,real(ϕ_pml[idxs]),m=:star,label="pml")
