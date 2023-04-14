using HarmonicWaterWaves
using LinearAlgebra
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW

# gr()

a, l = 10, 10
h    = 0.1
q    = 3
θ    = π/4
d    = 5

p    = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.set_depth!(tank,d)

HarmonicWaterWaves.add_freesurface!(tank,0,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

HarmonicWaterWaves.add_bottom!(tank,0,a)
HarmonicWaterWaves.add_bottom!(tank,a,a+l)

s1 = WPB.line(WW.Point2D(0,0),WW.Point2D(0,-d))
HarmonicWaterWaves.add_obstacles!(tank,WPB.Domain(s1))

HarmonicWaterWaves.add_orthogonal_pml!(tank;a,θ)

HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

plot(tank)

HarmonicWaterWaves.assemble_operators!(tank)

ϕ, dϕ = HarmonicWaterWaves.plane_wave(tank)
ψ, dψ = HarmonicWaterWaves.evanescent_wave(tank;s=1)
ϕi = (dof) -> 0.1*ϕ(dof) + ψ(dof)
dϕi = (dof) -> 0.1*dϕ(dof) + dψ(dof)

τ = WW.pml(tank)

k = WW.impedance(tank)
quad = tank.quad

idxs_free   = WPB.dom2qtags(quad,WW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(quad,WW.bottom(tank))
idxs_wall   = WPB.dom2qtags(quad,WW.obstacles(tank))

# ee1  = [dϕi(dof) - k*ϕi(dof) for dof in tank.quad.qnodes[idxs_free]]
# ee2  = [dϕi(dof) for dof in tank.quad.qnodes[idxs_bottom]]
# @show norm(ee1,Inf), norm(ee2,Inf)
##

function f(dof)
    x = dof.coords
    if x[1] == 0
        dϕi(dof)
    else
        return zero(ComplexF64)
    end
end

ϕ_pml = HarmonicWaterWaves.solve(tank,f)
ϕ_exact = [ϕi(dof) for dof in quad.qnodes]

x̃ = [τ(dof) for dof in quad.qnodes]
ϕ_pml_exact = [ϕi(x) for x in x̃]

#=
    Plot exact solution vs. PML solution on the free surface
=#
##
scale_size = 1.0
Plots.scalefontsizes(scale_size)
Γf = WW.freesurface(tank)
xx = [dof.coords[1] for dof in quad.qnodes]
# fig = plot(tank,aspect_ratio=4)
fig = plot(;legendfontsize=12)
labeled = false
for Γ in Γf
    I = WPB.dom2qtags(quad,Γ)
    plot!(fig,xx[I],real(ϕ_exact[I]),label= labeled ? "" : "exact",color=:blue,lw=2,ls=:dash)
    plot!(fig,xx[I],real(ϕ_pml_exact[I]),label= labeled ? "" : "exact pml",color=:green,lw=2,ls=:dash)
    plot!(fig,xx[I],real(ϕ_pml.vals[I]),label= labeled ? "" : "approximate pml",color=:red,ls=:solid)
    labeled = true
end
display(fig)
savefig("surface_wave_two_modes.png")
Plots.scalefontsizes()
