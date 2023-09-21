using LinearAlgebra
using CairoMakie
using LaTeXStrings
import WavePropBase as WPB
import HarmonicWaterWaves as HWW
using Unitful

HWW.set_makie_theme!()

impedances = [4 8 16 32]
ll   = 0.01:0.5:5.01 # PML size in wavelengths
dict = Dict()
pml_strength = 1
scales = Float64[]

for ν in impedances
    tank = HWW.WaveTank(; frequency = sqrt(ν)*u"Hz", depth = 1*u"m", gravity = 1*u"m/s^2")
    λ    = HWW.wavelength(tank)
    k    = HWW.wavenumber(tank)
    γ    = HWW.evanescent_wavenumber(tank)
    @assert HWW.impedance(tank) ≈ ν
    pml_start = 2*λ
    meshsize  = λ/10
    er = Float64[]
    for l in ll
        domain_end = pml_start + l*λ
        stretch_start = (γ*ν*domain_end + k*pml_start) / (k + γ*ν)
        @assert stretch_start > pml_start
        @assert stretch_start < domain_end
        τ = HWW.OrthogonalPML(; pml_start, stretch_start, pml_strength, stretch_strength=ν)
        empty!(tank)
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
    xlabel=L"(M-a) / \lambda",
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

l1_factor = π/(π+2*pml_strength)
ref = exp.(-2π.*l1_factor*ll) # reference slope is exp(-k*x), but we normalized k by λ = 2π/k
iref = findlast(l->l<3,ll)
lines!(ax,ll, dict[last(impedances)].er[iref]/ref[iref]*ref,
    label=L"\mathcal{O}(e^{-k \pi/(\pi+2) M})",
    linestyle=:dash,
    color=Cycled(length(impedances))
)

# ref = exp.(-2π/5 .*ll) # reference slope is exp(-k*x), but we normalized k by λ = 2π/k
# iref = findlast(l->l<3,ll)
# lines!(ax,ll, dict[last(Froud_sq)].er[iref]/ref[iref]*ref,
#     label=L"\mathcal{O}(e^{-2 \pi/2 M/\lambda})",
#     linestyle=:dash,
# )

leg = axislegend(ax;
    position=:lb
)
save("paper/figures/convergence_pml_stretch_planewave_vary_depth.pdf",fig)
