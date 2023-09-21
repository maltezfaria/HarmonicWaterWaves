using LinearAlgebra
using CairoMakie
using LaTeXStrings
import WavePropBase as WPB
import HarmonicWaterWaves as HWW
using Unitful

HWW.set_makie_theme!()

#
tank = HWW.WaveTank(; frequency = sqrt(9.81)u"Hz", depth = 1u"m", gravity = 9.81u"m/s^2")
L = HWW.lengthscale(tank)
λ = HWW.wavelength(tank)
k = HWW.wavenumber(tank)
γ₁ = HWW.evanescent_wavenumber(tank)

@info "k  = $k"
@info "γ₁ = $(γ₁)"

# construct the PML function
pml_start = 2 * λ
domain_end = pml_start + λ
τ = HWW.OrthogonalPML(; pml_start, pml_strength = 1)
HWW.add_pml!(tank, τ)

# add the wavemaker domain
HWW.add_wavemaker!(tank, domain_end)

# discretize domain
meshsize = λ / 10
HWW.discretize!(tank; meshsize, order = :high)
##

# assemble linear operators
HWW.assemble_operators!(tank)

# exact solution
ϕ, dϕ = HWW.plane_wave(tank)
ψ, dψ = HWW.evanescent_wave(tank; s = 1)
ϕi = (dof) -> ϕ(dof) + ψ(dof)
dϕi = (dof) -> dϕ(dof) + dψ(dof)

# boundary condition
ν = HWW.impedance(tank)
quad = tank.quad
idxs_free = WPB.dom2qtags(quad, HWW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(quad, HWW.bottom(tank))
idxs_wall = WPB.dom2qtags(quad, HWW.obstacles(tank))

f = zeros(ComplexF64, length(tank.quad.qnodes))
for i in idxs_wall ∪ idxs_bottom
    f[i] = dϕi(tank.quad.qnodes[i])
end
@info norm([dϕi(q) + ν * ϕi(q) for q in tank.quad.qnodes[idxs_free]], Inf)
@info norm([dϕi(q) for q in tank.quad.qnodes[idxs_bottom]], Inf)

# solution
ϕ_pml = HWW.solve!(tank, f)
## Plot exact solution and its analytic extension
Γf = HWW.freesurface(tank)
xx = [dof.coords[1] for dof in quad.qnodes]
# fig = plot(tank,aspect_ratio=4)
fig = Figure()
ax = Axis(fig[1, 1]; xlabel = L"x_1/\lambda", limits = (0, domain_end / λ, -1.1, 2.1))
xrange = 0:0.01:domain_end
ϕ_exact = [ϕi((x, 0)) for x in xrange]
ϕ_pml_exact = [ϕi(τ((x, 0))) for x in xrange]
lines!(
    ax,
    xrange ./ λ,
    real.(ϕ_exact);
    label = L"$\mathrm{Re}(\varphi)$",
    color = Cycled(1),
)
lines!(
    ax,
    xrange ./ λ,
    real(ϕ_pml_exact);
    label = L"$\mathrm{Re}(\tilde{\varphi})$",
    color = Cycled(2),
)

# plot the numerical solution
for Γ in Γf
    idxs = WPB.dom2qtags(quad, Γ)
    # lines!(ax, xx[I], real(ϕ_exact[I]), label=L"$\mathrm{Re}(\varphi)$", color=Cycled(1))
    scatter!(
        ax,
        xx[idxs] / λ,
        real(ϕ_pml.vals[idxs]);
        label = L"$\mathrm{Re}(\tilde{\varphi}_{M,h})$",
        color = Cycled(3),
    )
end

vlines!(ax, [pml_start] / λ; linestyle = :dash, color = (:black, 0.5))

leg = axislegend(ax; unique = true, merge = true, position = :rt)
save("paper/figures/wavemaker_modal_solution.pdf", fig)
