using HarmonicWaterWaves
using LinearAlgebra
using Nystrom
using Nystrom: ParametricSurfaces
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW

gr()

a, l = 10, 2
h    = 0.1
q    = 3
θ    = π/4
d    = 1

p    = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.set_depth!(tank,d)

HarmonicWaterWaves.add_freesurface!(tank,0,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

HarmonicWaterWaves.add_bottom!(tank,0,a)
HarmonicWaterWaves.add_bottom!(tank,a,a+l)

s1 = ParametricSurfaces.line(WW.Point2D(0,0),WW.Point2D(0,-d))
HarmonicWaterWaves.add_obstacles!(tank,Domain(s1))

HarmonicWaterWaves.add_pml!(tank;a,θ)

plot(tank)

HarmonicWaterWaves.discretize!(tank;meshsize=h,order=q)

HarmonicWaterWaves.assemble_operators!(tank)

ϕi, dϕi = HarmonicWaterWaves.plane_wave(tank)


τ = WW.pml(tank)

k = WW.impedance(tank)
quad = tank.quad

idxs_free   = Nystrom.dom2dof(quad,WW.freesurface(tank))
idxs_bottom = Nystrom.dom2dof(quad,WW.bottom(tank))
idxs_wall   = Nystrom.dom2dof(quad,WW.obstacles(tank))

# ff  = [ϕi(dof) - k*dϕi(dof) for dof in tank.quad.dofs[idxs_free]]
# ff  = [dϕi(dof) for dof in tank.quad.dofs[idxs_bottom]]

function f(dof)
    x = dof.coords
    if x[1] == 0
        dϕi(dof)
    else
        return zero(ComplexF64)
    end
end

ϕ_pml = HarmonicWaterWaves.solve(tank,f)
ϕ_exact = [ϕi(dof) for dof in quad.dofs]

x̃ = [τ(dof) for dof in quad.dofs]
ϕ_pml_exact = [ϕi(x) for x in x̃]

#=
    Plot exact solution vs. PML solution on the free surface
=#
##
scale_size = 1.5
Plots.scalefontsizes(scale_size)
Γf = WW.freesurface(tank)
xx = [dof.coords[1] for dof in quad.dofs]
# fig = plot(tank,aspect_ratio=4)
fig = plot(;legendfontsize=12)
labeled = false
for Γ in Γf
    I = Nystrom.dom2dof(quad,Γ)
    plot!(fig,xx[I],real(ϕ_exact[I]),label= labeled ? "" : "exact",color=:blue,lw=2,ls=:dash)
    plot!(fig,xx[I],real(ϕ_pml_exact[I]),label= labeled ? "" : "exact pml",color=:green,lw=2,ls=:dash)
    plot!(fig,xx[I],real(ϕ_pml.vals[I]),label= labeled ? "" : "approximate pml",color=:red,ls=:solid)
    labeled = true
end
display(fig)
savefig("surface_wave.png")
Plots.scalefontsizes()
