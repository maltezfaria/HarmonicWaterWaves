using Test
using LinearAlgebra
using CairoMakie
using LaTeXStrings
import WavePropBase as WPB
import HarmonicWaterWaves as HWW
using Unitful

HWW.set_makie_theme!()

ν     = 4
depth = 1 * u"m"
tank  = HWW.WaveTank(; depth, frequency = sqrt(ν) * u"Hz", gravity = 1 * u"m/s^2")

λ = HWW.wavelength(tank)

pml_start = 3λ
pml_length = 2λ
domain_end = pml_start + pml_length

# add pml
τ = HWW.OrthogonalPML(; pml_start)
HWW.add_pml!(tank, τ)

# piercing obstacles
r = 0.25 # radius of obstacles
δ = 2 * r + 1 # distance between center of obstacles
x₁ = (-δ / 2, 0.0) # center of first obstacle
x₂ = (δ / 2, 0.0) # center of first obstacle
circ1 = let r = r, x₀ = x₁
    WPB.ParametricEntity(0, π) do (u,)
        return (-r * cos(u), -r * sin(u)) .+ x₀
    end
end
circ2 = let r = r, x₀ = x₂
    WPB.ParametricEntity(0, π) do (u,)
        return (-r * cos(u), -r * sin(u)) .+ x₀
    end
end
obs = WPB.Domain([circ1, circ2])
HWW.add_obstacles!(tank, obs)

# free surface. Care is taken to exclude the piercing the obstacles, and to have
# a grid poin on the PML boundary
HWW.add_freesurface!(tank, -domain_end, -pml_start)
HWW.add_freesurface!(tank, -pml_start, -δ / 2 - r)
HWW.add_freesurface!(tank, -δ / 2 + r, δ / 2 - r)
HWW.add_freesurface!(tank, δ / 2 + r, pml_start)
HWW.add_freesurface!(tank, pml_start, domain_end)

# bottom
HWW.add_bottom!(tank, -domain_end, -pml_start)
HWW.add_bottom!(tank, -pml_start, pml_start)
HWW.add_bottom!(tank, pml_start, domain_end)

ϕi, dϕi = HWW.plane_wave(tank)

P     = 5
ppws  = [4 * 2^i for i in 0:6]
sols  = Dict()
Xtest = [(x, y) for x in 1.5*λ:0.1:2.5λ, y in -0.75:0.1:-0.25]

for ppw in ppws
    # mesh the obstacles with a smaller mesh size due to their curvature
    h = P * λ / ppw
    meshsize = h
    # create mesh
    HWW.discretize!(tank; meshsize, order = P)
    idxs_free   = WPB.dom2qtags(tank.quad, HWW.freesurface(tank))
    idxs_bottom = WPB.dom2qtags(tank.quad, HWW.bottom(tank))
    idxs_obs    = WPB.dom2qtags(tank.quad, HWW.obstacles(tank))
    f           = zeros(ComplexF64, length(tank.quad.qnodes))
    for i in idxs_obs ∪ idxs_bottom
        f[i] = -dϕi(tank.quad.qnodes[i])
    end
    # sanity checks
    ee1 = [dϕi(dof) + HWW.impedance(tank) * ϕi(dof) for dof in tank.quad.qnodes[idxs_free]]
    if norm(ee1, Inf) > 1e-10
        @warn "Incident wave does not satify surface bottom conditions" norm(ee1, Inf),
        norm(ee2, Inf)
    end

    # solve
    HWW.assemble_operators!(tank)
    ϕ_pml = HWW.solve!(tank, f)

    # compute solution on test domain
    sol = HWW.solution(tank)
    sols[ppw] = sol.(Xtest)
end

# now compute the self-convergence errors
sol_ref = sols[ppws[end]]
ee = Float64[]
for ppw in ppws[1:end-1] # skipe last, which is the reference solution
    sol = sols[ppw]
    err = norm(sol - sol_ref, Inf) / norm(sol_ref, Inf)
    push!(ee, err)
end

## plot error
fig = Figure()
ax  = Axis(fig[1, 1]; xscale = log10, yscale = log10, xticks = ppws, xminorticksvisible = false, yminorticksvisible = false, xlabel = "points per wavelength", ylabel = "error")

scatterlines!(ax, ppws[1:end-1], ee; label = L"P = %$P")
ref = 1 ./ ppws[1:end-1] .^ (big(P))
iref = length(ppws) - 1
lines!(ax, ppws[1:end-1], ee[iref] / ref[iref] * ref; label = L"\mathcal{O}(h^{-%$P})", linestyle = :dash)
# ylims!(ax,1e-12,1)
xlims!(ax, first(ppws), ppws[end-1])
axislegend(ax; position = :rt)

save("paper/figures/mesh_convergence_double_piercing.pdf", fig)
