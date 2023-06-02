using HarmonicWaterWaves
using LinearAlgebra
using Plots
using LaTeXStrings
import WavePropBase as WPB
import HarmonicWaterWaves as WW

gr()

savefig_ = false

pml_start     = 10
pml_length    = 10
domain_end    = stretch_start + stretch_length
meshsize      = 0.5
order         = 7 # method order
qorder        = 2*(order-1)-1 # quadrature order
depth         = 1

p    = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.set_depth!(tank,depth)

HarmonicWaterWaves.add_freesurface!(tank,0,pml_start)
HarmonicWaterWaves.add_freesurface!(tank,pml_start,domain_end)

HarmonicWaterWaves.add_bottom!(tank,0,pml_start)
HarmonicWaterWaves.add_bottom!(tank,pml_start,domain_end)

side = WPB.line(WW.Point2D(0,0),WW.Point2D(0,-depth))

HarmonicWaterWaves.add_obstacles!(tank,WPB.Domain(side))

pml = WW.OrthogonalPML(;a=pml_start)

HarmonicWaterWaves.add_pml!(tank,pml)
HarmonicWaterWaves.discretize!(tank;meshsize,qorder)

HarmonicWaterWaves.assemble_operators!(tank)

ϕ, dϕ = HarmonicWaterWaves.plane_wave(tank)
ψ, dψ = HarmonicWaterWaves.evanescent_wave(tank;s=1)
ϕi = (dof)  -> ϕ(dof) + ψ(dof)
dϕi = (dof) -> dϕ(dof) + dψ(dof)

τ = WW.pml(tank)

k = WW.impedance(tank)
quad = tank.quad

idxs_free   = WPB.dom2qtags(quad,WW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(quad,WW.bottom(tank))
idxs_wall   = WPB.dom2qtags(quad,WW.obstacles(tank))

function f(dof)
    x = dof.coords
    if x[1] == 0
        dϕi(dof)
    else
        return zero(ComplexF64)
    end
end

ϕ_pml = HarmonicWaterWaves.solve!(tank,f)
ϕ_exact = [ϕi(dof) for dof in quad.qnodes]

x̃ = [τ(dof) for dof in quad.qnodes]
ϕ_pml_exact = [ϕi(x) for x in x̃]

#=
    Plot exact solution vs. PML solution on the free surface
=#
##
Γf = WW.freesurface(tank)
xx = [dof.coords[1] for dof in quad.qnodes]
# fig = plot(tank,aspect_ratio=4)
fig = plot(;legendfontsize=12,xtickfontsize=10,ytickfontsize=10,xlabel=L"x")
labeled = false
for Γ in Γf
    I = WPB.dom2qtags(quad,Γ)
    plot!(fig,xx[I],real(ϕ_exact[I]),label= labeled ? "" : L"$\Re(\varphi)$",color=:blue,lw=2,ls=:dot)
    plot!(fig,xx[I],real(ϕ_pml_exact[I]),label= labeled ? "" : L"$\Re(\tilde{\varphi})$",color=:green,lw=2,ls=:solid)
    plot!(fig,xx[I],real(ϕ_pml.vals[I]),label= labeled ? "" : L"$\Re(\tilde{\varphi}_{\ell,h})$",color=:red,ls=:dash,lw=2)
    labeled = true
end

display(fig)

savefig_ && savefig(fig,"paper/figures/wavemaker_modal_solution.pdf")


##
# sol = WW.solution(tank)
# xx = 0:0.025:a+l
# yy = -d:0.025:0
# u = [sol((x,y)) for y in yy, x in xx]

# heatmap(xx,yy,real.(u))

# x0  = WPB.Point2D(3,-d/3)
# sol(x0) + ϕi(x0)
