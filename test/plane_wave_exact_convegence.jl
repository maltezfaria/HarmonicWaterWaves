using HarmonicWaterWaves
using LinearAlgebra
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW
using Roots
using LaTeXStrings

gr()
h    = 0.1
a    = 1
q    = 6
θ    = π/4
d    = 20

ll = 0.1:0.5:10 |> collect
er = []
p = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1,depth=d)
for l in ll
    tank = HarmonicWaterWaves.WaveTank(parameters=p)
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
    ϕi = (dof) -> 0.1*ϕ(dof)   + ψ(dof)
    dϕi = (dof) -> 0.1*dϕ(dof) + dψ(dof)

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

    ϕ_pml = HarmonicWaterWaves.solve(tank,f)
    ϕ_exact = [ϕi(dof) for dof in quad.qnodes]

    x̃ = [τ(dof) for dof in quad.qnodes]
    ϕ_pml_exact = [ϕi(x) for x in x̃]

    idxs_test = filter(i -> 0.2*a < quad.qnodes[i].coords[1] < 0.8*a,idxs_free)
    push!(er,norm(ϕ_exact[idxs_test] - ϕ_pml[idxs_test],Inf))
    # @info l,er[end]
end

##
fig = plot(yscale=:log10,xlabel=L"\ell",ylabel="error",m=:x,ylims=(5e-3,1e1))
plot!(fig,ll,er,yscale=:log10,m=:x,label="d=$d")
# plot reference exponential decay
β = sin(θ)*WW.wavenumber(p)
γ = WW.evanescent_wavenumber(p)
ref = exp.(-β.*ll)./ll
iref = findlast(l->l<15,ll)
plot!(fig,ll,er[iref]/ref[iref]*ref,label=L"\frac{e^{-k \sin(\theta)\ell}}{\ell}",ls=:dash)
ref = exp.(-β.*ll)
plot!(fig,ll,er[iref]/ref[iref]*ref,label=L"e^{-k \sin(\theta)\ell}",ls=:dash)

@show β, γ

fig



# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")
# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")

# fig2 = deepcopy(fig)

savefig(fig,"convergence_pml.png")