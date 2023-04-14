using HarmonicWaterWaves
using LinearAlgebra
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW
using Roots
using LaTeXStrings

gr()
h    = 0.1
a    = 10
q    = 8
θ    = π/4
d = 0.25

fig = plot(yscale=:log10,xlabel=L"\ell",ylabel="error",m=:x)

ll = 0.1:0.5:8 |> collect
er = []
γ = find_zero(x->x*tan(x*d)-1,(0,π/(2d))) * cos(θ)
p = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1,depth=d)
β = sin(θ)*WW.wavenumber(p)
@info γ,β
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

    ϕi, dϕi = HarmonicWaterWaves.plane_wave(tank)

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
            dϕi(dof)
        else
            return zero(ComplexF64)
        end
    end

    ϕ_pml = HarmonicWaterWaves.solve(tank,f)
    ϕ_exact = [ϕi(dof) for dof in quad.qnodes]

    x̃ = [τ(dof) for dof in quad.qnodes]
    ϕ_pml_exact = [ϕi(x) for x in x̃]

    idxs_test = filter(i -> 0.0*a < quad.qnodes[i].coords[1] < 1.0*a,idxs_free)
    push!(er,norm(ϕ_exact[idxs_test] - ϕ_pml[idxs_test],Inf))
    # @info l,er[end]
end
plot!(fig,ll,er,yscale=:log10,m=:x,label="d=$d")
ref = exp.(-β.*ll)
plot!(fig,ll,er[end]/ref[end]*ref,label=L"e^{-k \cos(\theta)\ell}",ls=:dash)

# fig

# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")
# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")

# fig2 = deepcopy(fig)

# savefig(fig,"convergence_pml.png")
