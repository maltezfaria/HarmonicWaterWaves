using HarmonicWaterWaves
using LinearAlgebra
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW
using Roots
using LaTeXStrings

gr()
a     = 1
order  = 7 # method order
n = order - 1 # number of nodes per element
q    = 2*n-1 # quadrature order
d    = 1

er_dict = Dict()
p = HarmonicWaterWaves.Parameters(frequency=sqrt(4),gravity=1,depth=d)
λ = WW.wavelength(p)
l = 2λ

ppw = [10,20,40]
cc = 0.1:0.2:3 |> collect

for ppw in ppw
    h = n * λ / ppw
    er = []
    for c in cc
        tank = HarmonicWaterWaves.WaveTank(parameters=p)
        HarmonicWaterWaves.add_freesurface!(tank,0,a)
        HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

        HarmonicWaterWaves.add_bottom!(tank,0,a)
        HarmonicWaterWaves.add_bottom!(tank,a,a+l)

        s1 = WPB.line(WW.Point2D(0,0),WW.Point2D(0,-d))
        HarmonicWaterWaves.add_obstacles!(tank,WPB.Domain(s1))

        pml = WW.OrthogonalPML(;a,c)
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
        ϕ_exact = [ϕi(dof) for dof in quad.qnodes]

        x̃ = [τ(dof) for dof in quad.qnodes]
        ϕ_pml_exact = [ϕi(x) for x in x̃]

        idxs_test = filter(i -> 0.0*a < quad.qnodes[i].coords[1] < 1.0*a,idxs_free)
        push!(er,norm(ϕ_exact[idxs_test] - ϕ_pml[idxs_test],Inf))
    end
    er_dict[ppw] = er
end

##
default(legendfontsize=12,xlabelfontsize=12,ylabelfontsize=12,guidefontsize=12,lw=2)
fig = plot(xlabel="pml attenuation",m=:x,ylabel="error")

for ppw in ppw
    plot!(fig, cc, er_dict[ppw], yscale=:log10, m=:x, label=L"ppw=%$ppw")
end

savefig(fig,"paper/figures/pml_attenuation.pdf")

fig

# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")
# plot(ll,er,yscale=:log10,xlabel="L",ylabel="error",m=:x,label="h=$h")

# fig2 = deepcopy(fig)

# savefig(fig,"convergence_pml.png")
