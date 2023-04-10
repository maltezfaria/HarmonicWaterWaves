using HarmonicWaterWaves
using LinearAlgebra
using Nystrom
using Nystrom: ParametricSurfaces
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW

gr()

ll = 0.1:0.1:8 |> collect
er = []

for l in ll
    a    = 10
    h    = 0.1
    q    = 3
    θ    = π/4
    d    = 0.25

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

    idxs_test = filter(i -> 0.0*a < quad.dofs[i].coords[1] < 1.0*a,idxs_free)
    push!(er,norm(ϕ_exact[idxs_test] - ϕ_pml[idxs_test],Inf))
    @info l,er[end]
end

plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=0.1")

# fig2 = deepcopy(fig)

# savefig(fig,"convergence_pml.png")
