using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 4,4π
h = 0.1
q = 3
θ = π/4

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

k = WW.impedance(tank)
op  = WW.InfiniteDepthWaterWaves(;k)
G   = WPB.SingleLayerKernel(op)
dG  = WPB.DoubleLayerKernel(op)

ϕe   = source -> G(xₛ,source)
dϕe  = source -> dG(xₛ,source)



## Solve with proper Green function
t_green = @elapsed begin
    Γs         = WW.obstacles(tank)
    quad_obs   = WPB.NystromMesh(view(tank.mesh,Γs);qorder=q)
    op         = WW.InfiniteDepthWaterWaves(;k)
    Sop        = WPB.SingleLayerOperator(op,quad_obs)
    Dop        = WPB.DoubleLayerOperator(op,quad_obs)
    S₀         = Sop |> Matrix
    D₀         = Dop |> Matrix
    δS         = WPB.hcubature_correction(Sop;max_dist=5*h,maxevals=20)
    δD         = WPB.hcubature_correction(Dop;max_dist=5*h,maxevals=20)
    S = S₀ + δS
    D = D₀ + δD
    f          = [dϕe(dof) for dof in quad_obs.qnodes]
    ϕ_green = (I/2 + D)\(S*f)
end

xx   = [atan(dof.coords[2]-x₀[2],dof.coords[1]-x₀[1]) for dof in quad_obs.qnodes]
yy   = [ϕe(dof) for dof in quad_obs.qnodes]

@info t_green, length(quad_obs.qnodes), norm(yy-ϕ_green)

scatter(xx,real(yy),label="exact",m=:cross,ms=10)
scatter!(xx,real(ϕ_green),label="green",m=:circle)
##

# Solve with PML Green function
# plot(tank)

t_pml = @elapsed begin
    quad = tank.quad
    dofs = quad.qnodes
    HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
    f = [dof.coords[2] == 0 ? dϕe(dof) - k*ϕe(dof) : dϕe(dof) for dof in dofs]
    ϕ_pml = HarmonicWaterWaves.solve(tank,f)
end

idxs = WPB.dom2qtags(quad,Γs)

@info t_pml, length(quad.qnodes), norm(yy-ϕ_pml[idxs])

scatter!(xx,real(ϕ_pml[idxs]),m=:star,label="pml")
