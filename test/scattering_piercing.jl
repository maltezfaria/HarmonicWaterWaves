using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 2,2π
h = 0.05
q = 3
θ = π/4
d = 2

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

# create piercing obstacle
r = 0.5 # radius of obstacle
x₀  = SVector(0,0.0) # center of obstacle
obs = let r = r, x₀ = x₀
    WPB.ParametricEntity(0,π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₀
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

# free surface
HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,-r)
# HarmonicWaterWaves.add_freesurface!(tank,-r,r)
HarmonicWaterWaves.add_freesurface!(tank,r,a)
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

# boundary conditions
ϕi, dϕi = HarmonicWaterWaves.plane_wave(tank)

idxs_free   = WPB.dom2qtags(tank.quad,WW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(tank.quad,WW.bottom(tank))
idxs_obs   = WPB.dom2qtags(tank.quad,WW.obstacles(tank))

f = zeros(ComplexF64, length(tank.quad.qnodes))
for i in idxs_obs
    f[i] = -dϕi(tank.quad.qnodes[i])
end

# sanity checks
ee1  = [dϕi(dof) - WW.impedance(tank)*ϕi(dof) for dof in tank.quad.qnodes[idxs_free]]
ee2  = [dϕi(dof) for dof in tank.quad.qnodes[idxs_bottom]]
if norm(ee1,Inf) > 1e-10 || norm(ee2,Inf) > 1e-10
    @warn "Incident wave does not satify surface and/or bottom conditions" norm(ee1,Inf), norm(ee2,Inf)
end

# solve
HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
ϕ_pml = HarmonicWaterWaves.solve(tank,f)

##
scale_size = 1.0
Plots.scalefontsizes(scale_size)
Γf = WW.freesurface(tank)
xx = [dof.coords[1] for dof in tank.quad.qnodes]
# fig = plot(tank,aspect_ratio=4)
fig = plot(;legendfontsize=12)
labeled = false
for Γ in Γf
    I = WPB.dom2qtags(tank.quad,Γ)
    plot!(fig,xx[I],real(ϕi.(tank.quad.qnodes[I]) +  ϕ_pml.vals[I]),label= labeled ? "" : "approximate pml",color=:red,ls=:solid)
    labeled = true
end
plot!(fig,tank)
display(fig)
