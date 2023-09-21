using LinearAlgebra
using CairoMakie
using LaTeXStrings
import HarmonicWaterWaves as HWW
import WavePropBase as WPB
using Unitful

HWW.set_makie_theme!()

impedances = [1,2,4,8]
ll   = 0.5:0.5:5 # PML size in wavelengths
dict = Dict()

for ν in impedances
    tank = HWW.WaveTank(; frequency = sqrt(ν)*u"Hz", depth = 1*u"m", gravity = 1*u"m/s^2")
    L    = HWW.lengthscale(tank)
    λ    = HWW.wavelength(tank)
    k    = HWW.wavenumber(tank)
    γ    = HWW.evanescent_wavenumber(tank)
    # dimensionless number
    @assert HWW.impedance(tank) ≈ ν
    pml_start = 2*λ
    τ = HWW.OrthogonalPML(; pml_start)
    meshsize  = λ/10
    er = Float64[]
    for l in ll
        empty!(tank)
        domain_end = pml_start + l*λ
        HWW.add_pml!(tank,τ)
        HWW.add_wavemaker!(tank,domain_end)
        HWW.discretize!(tank;meshsize,order=:high)
        HWW.assemble_operators!(tank)
        ϕ, dϕ = HWW.plane_wave(tank)
        ψ, dψ = HWW.evanescent_wave(tank;s=1)
        ϕi = (dof)  -> ϕ(dof) + ψ(dof)
        dϕi = (dof) -> dϕ(dof) + dψ(dof)
        quad = tank.quad
        idxs_free   = WPB.dom2qtags(quad,HWW.freesurface(tank))
        idxs_bottom = WPB.dom2qtags(quad,HWW.bottom(tank))
        idxs_wall   = WPB.dom2qtags(quad,HWW.obstacles(tank))
        f = zeros(ComplexF64, length(tank.quad.qnodes))
        for i in idxs_wall
            f[i] = dϕi(tank.quad.qnodes[i])
        end
        ϕ_pml = HWW.solve!(tank,f)
        ϕ_exact = [ϕi(τ(dof)) for dof in quad.qnodes]
        x̃ = [τ(dof) for dof in quad.qnodes]
        ϕ_pml_exact = [ϕi(x) for x in x̃]
        push!(er,norm(ϕ_exact[idxs_free] - ϕ_pml[idxs_free],Inf))
        dict[ν] = (er=er,γ=γ,λ=λ,k=k)
    end
end

## Plot the error for various depths
fig = Figure()
ax  = Axis(fig[1,1],
    yscale=log10,
    xlabel=L"(M-a)/\lambda",
    ylabel=L"E_{M,h}",
    limits=(0,last(ll),1e-12,1),
    yminorticksvisible=false,
)

for ν in impedances
    vals = dict[ν]
    sc = scatterlines!(ax,ll, vals.er,
        marker=:circle,
        label=L"\nu = %$(ν)"
    )
end

ref = exp.(-2π.*ll) # reference slope is exp(-k*x), but we normalized k by λ = 2π/k
iref = findlast(l->l<3,ll)
lines!(ax,ll, dict[1].er[iref]/ref[iref]*ref,
    marker=:circle,
    label=L"\mathcal{O}(e^{-M})",
    linestyle=:dash,
)

# also add the evanescent wave reference slope for the deepest case
γ = dict[impedances[end]].γ
λ = 2π / dict[impedances[end]].k
ref = exp.(-γ*λ.*ll) # exp(-γ*M) since M = λ*ll
iref = findlast(l->l<3,ll)
lines!(ax,ll, dict[impedances[end]].er[iref]/ref[iref]*ref,
    color=Cycled(length(impedances)),
    label=L"\mathcal{O}(e^{-\gamma_1 M})",
    linestyle=:dash,
)

leg = axislegend(ax;
    position=:lb
)
save("paper/figures/convergence_pml_planewave_vary_depth.pdf",fig)
