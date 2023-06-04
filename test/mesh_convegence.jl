using HarmonicWaterWaves
using LinearAlgebra
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW
using Roots
using LaTeXStrings

gr()
a      = 1
order  = 6 # method order
n      = order # number of nodes per element
q      = 2*n-1 # quadrature order
d      = 1

er = []
p = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1,depth=d)
λ = WW.wavelength(p)
l = 8λ

hh =[λ*1/2^i for i in 0:6]

for h in hh
    tank = HarmonicWaterWaves.WaveTank(parameters=p)
    HarmonicWaterWaves.add_freesurface!(tank,0,a)
    HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

    HarmonicWaterWaves.add_bottom!(tank,0,a)
    HarmonicWaterWaves.add_bottom!(tank,a,a+l)

    s1 = WPB.line(WW.Point2D(0,0),WW.Point2D(0,-d))
    HarmonicWaterWaves.add_obstacles!(tank,WPB.Domain(s1))

    pml = WW.OrthogonalPML(;a,c=1.0)
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

    function f(dof)
        x = dof.coords
        if x[1] == 0
            ComplexF64(dϕi(dof))
        else
            return zero(ComplexF64)
        end
    end

    ϕ_pml = HarmonicWaterWaves.solve!(tank,f)
    ϕ_exact = [ϕi(τ(dof)) for dof in quad.qnodes]

    x̃ = [τ(dof) for dof in quad.qnodes]
    ϕ_pml_exact = [ϕi(x) for x in x̃]

    # idxs_test = filter(i -> 0.0*a < quad.qnodes[i].coords[1] < 1.0*a,idxs_free)
    # push!(er,norm(ϕ_exact[idxs_test] - ϕ_pml[idxs_test],Inf))
    push!(er,norm(ϕ_exact - ϕ_pml,Inf))
    # @info h,er[end]
end

##
ppw = n * λ ./ hh
default(legendfontsize=12,xlabelfontsize=12,ylabelfontsize=12,guidefontsize=12,lw=2)
fig = plot(yscale=:log10,xscale=:log10,xlabel="points per wavelength",m=:x)
plot!(fig, ppw, er,yscale=:log10,m=:x,label=L"|| \varphi_{\ell,h} - \varphi ||")
# plot reference decay
# ref = 1 ./ ppw.^(order) .* (log.(ppw))
ref = 1 ./ ppw.^(order)
# ref = hh.^(order)
# iref = findlast(h->h<5,ll_norm)
iref = length(ppw) - 1
plot!(fig,ppw,er[iref]/ref[iref]*ref,label=L"\mathcal{O}(h^%$order \log(1/h)) ",ls=:dash)

xlabels = [string(Int(p)) for p in ppw]
xticks!(fig,(ppw,xlabels))

# # plot reference exponential decay
# if plot_evanescent
#     ref = exp.(-γ.*ll)
#     plot!(fig,ll_norm,er[iref]/ref[iref]*ref,label=L"e^{-\gamma_1 \ell}",ls=:dash)
# end
# iref = length(ll)

savefig(fig,"paper/figures/modal_convergence_meshsize_depth_$(d).pdf")

fig

# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")
# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")

# fig2 = deepcopy(fig)

# savefig(fig,"convergence_pml.png")
