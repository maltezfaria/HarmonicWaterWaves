using HarmonicWaterWaves
using LinearAlgebra
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW
using Roots
using LaTeXStrings

gr()
a    = 1
order = 7 # method order
q    = 2*(order-1)-1 # quadrature order
d    = 4

plot_evanescent = d != 1

er = []
p = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1,depth=d)
λ = WW.wavelength(p)

h    = λ / 30
b    = λ

ll = 0.1:1:4*λ |> collect

for l in ll

    tank = HarmonicWaterWaves.WaveTank(parameters=p)
    HarmonicWaterWaves.add_freesurface!(tank,0,a)
    HarmonicWaterWaves.add_freesurface!(tank,a,a+l)
    HarmonicWaterWaves.add_freesurface!(tank,a+l,a+l+b)

    HarmonicWaterWaves.add_bottom!(tank,0,a)
    HarmonicWaterWaves.add_bottom!(tank,a,a+l)
    HarmonicWaterWaves.add_bottom!(tank,a+l,a+l+b)

    s1 = WPB.line(WW.Point2D(0,0),WW.Point2D(0,-d))
    HarmonicWaterWaves.add_obstacles!(tank,WPB.Domain(s1))

    pml = WW.OrthogonalPML(;a=a,b=a+l,d=d)

    HarmonicWaterWaves.add_pml!(tank,pml)

    HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

    plot(tank)

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

    # ff  = [ϕi(dof) - k*dϕi(dof) for dof in tank.quad.dofs[idxs_free]]
    # ff  = [dϕi(dof) for dof in tank.quad.dofs[idxs_bottom]]

    function f(dof)
        x = dof.coords
        if x[1] == 0
            ComplexF64(dϕi(dof))
        else
            return zero(ComplexF64)
        end
    end

    ϕ_pml = HarmonicWaterWaves.solve!(tank,f)
    ϕ_exact = [ϕi(dof) for dof in quad.qnodes]

    x̃ = [τ(dof) for dof in quad.qnodes]
    ϕ_pml_exact = [ϕi(x) for x in x̃]

    idxs_test = filter(i -> 0.1*a < quad.qnodes[i].coords[1] < 0.8*a,idxs_free)
    push!(er,norm(ϕ_exact[idxs_test] - ϕ_pml[idxs_test],Inf))
    # @info l,er[end]
end

##
ll_norm = ll/λ
default(legendfontsize=12,xlabelfontsize=12,ylabelfontsize=12,guidefontsize=12)
fig = plot(yscale=:log10,xlabel=L"\ell / \lambda",m=:x,yticks=[10.0^(-i) for i in -1:12],ylims=(1e-13,0))
plot!(fig, ll_norm, er,yscale=:log10,m=:x,label=L"|| \varphi_{\ell,h} - \varphi ||_{\infty}")
# plot reference exponential decay
β = WW.wavenumber(p)
γ = WW.evanescent_wavenumber(p)
ref = exp.(-β.*ll)
iref = findlast(l->l<3,ll_norm)
plot!(fig,ll_norm,er[iref]/ref[iref]*ref,label=L"\mathcal{O}(e^{-k\ell})",ls=:dash)
# plot reference exponential decay
if plot_evanescent
    ref = exp.(-γ.*ll)
    plot!(fig,ll_norm,er[iref]/ref[iref]*ref,label=L"\mathcal{O}(e^{-\gamma_1 \ell})",ls=:dash)
end
# iref = length(ll)

savefig(fig,"paper/figures/convergence_pml_stretching_planewave_depth_$(d).pdf")

fig

# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")
# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")

# fig2 = deepcopy(fig)

# savefig(fig,"convergence_pml.png")
