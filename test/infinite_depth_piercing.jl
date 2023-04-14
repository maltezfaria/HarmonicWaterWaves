using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 1,4π
h = 0.4
q = 3
θ = π/4

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

r = 0.7
HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,-r)
# HarmonicWaterWaves.add_freesurface!(tank,-r,r)
HarmonicWaterWaves.add_freesurface!(tank,r,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

HarmonicWaterWaves.add_pml!(tank;a,θ)

x₀  = SVector(0,0.0) # center of obstacle
xₛ  = SVector(0.1,0.2) # center of source
obs = let r = r, x₀ = x₀
    WPB.ParametricEntity(0,π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₀
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

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

@info t_green, length(quad_obs.qnodes), norm(yy-ϕ_green,Inf)/norm(yy,Inf)

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

@info t_pml, length(quad.qnodes), norm(yy-ϕ_pml[idxs],Inf)/norm(yy,Inf)

scatter!(xx,real(ϕ_pml[idxs]),m=:star,label="pml")


scale_size = 1.0
Plots.scalefontsizes(scale_size)
Γf = WW.freesurface(tank)
xx = [dof.coords[1] for dof in quad.qnodes]
# fig = plot(tank,aspect_ratio=4)
fig = plot(;legendfontsize=12)
labeled = false
for Γ in Γf
    I = WPB.dom2qtags(quad,Γ)
    plot!(fig,xx[I],real(ϕ_pml.vals[I]),label= labeled ? "" : "approximate pml",color=:red,ls=:solid)
    labeled = true
end
display(fig)
