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

a, l = 2,4π
h = 0.1
q = 3
θ = π/4

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

HarmonicWaterWaves.add_pml!(tank;a,θ)

x₀  = SVector(1,-1.0) # center of obstacle
xₛ  = SVector(1.7,-1.1) # center of source
# xₛ  = x₀
obs = boundary(ParametricSurfaces.Disk(;radius=0.8,center=x₀,normal=:inside))
HarmonicWaterWaves.add_obstacles!(tank,obs)

HarmonicWaterWaves.discretize!(tank;meshsize=h,order=q)

k = WW.impedance(tank)
op  = WW.InfiniteDepthWaterWaves(;k)
G   = SingleLayerKernel(op)
dG  = DoubleLayerKernel(op)
ϕe   = source -> G(xₛ,source)
dϕe  = source -> dG(xₛ,source)

## Solve with proper Green function
t_green = @elapsed begin
    Γs         = WW.obstacles(tank)
    quad_obs   = NystromMesh(view(tank.mesh,Γs);order=q)
    op         = WW.InfiniteDepthWaterWaves(;k)
    S          = SingleLayerOperator(op,quad_obs) |> Nystrom.assemble_gk
    D          = DoubleLayerOperator(op,quad_obs) |> Nystrom.assemble_gk
    f          = [dϕe(dof) for dof in quad_obs.dofs]
    ϕ_green = (I/2 + D)\(S*f)
end

xx   = [atan(dof.coords[2]-x₀[2],dof.coords[1]-x₀[1]) for dof in quad_obs.dofs]
yy   = [ϕe(dof) for dof in quad_obs.dofs]

@info t_green, length(quad_obs.dofs), norm(yy-ϕ_green)

scatter(xx,real(yy),label="exact",m=:cross,ms=10)
scatter!(xx,real(ϕ_green),label="green",m=:circle)
##

## Solve with PML Green function
# plot(tank)

t_pml = @elapsed begin
    quad = tank.quad
    dofs = quad.dofs
    HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
    f = [dof.coords[2] == 0 ? dϕe(dof) - k*ϕe(dof) : dϕe(dof) for dof in dofs]
    ϕ_pml = HarmonicWaterWaves.solve(tank,f)
end

idxs = Nystrom.dom2dof(quad,Γs)

@info t_pml, length(quad.dofs), norm(yy-ϕ_pml[idxs])

scatter!(xx,real(ϕ_pml[idxs]),m=:star,label="pml")
