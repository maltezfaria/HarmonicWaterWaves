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
pml_length = λ

# add pml
τ = HWW.OrthogonalPML(; pml_start)
HWW.add_pml!(tank, τ)

# waveguide
depth_func = let λ = λ
    (s) -> -1 + 0.5 * exp(-s^2 / 4) * sin(2π * s / λ)
end
HWW.add_waveguide!(tank, pml_start + pml_length, depth_func)

# obstacles
jelly1 = HWW.jellyfish(1 / 4, (-0.5, -0.4), 0)
jelly2 = HWW.jellyfish(1 / 4, (-2, -0.5), π / 6)
jelly3 = HWW.jellyfish(1 / 4, (1, -0.5), -π / 6)
obs = WPB.Domain([jelly1, jelly2, jelly3])
obs = WPB.Domain([jelly1])
HWW.add_obstacles!(tank, obs)

ϕi, dϕi = HWW.plane_wave(tank)

P     = 5
ppws  = [4 * 2^i for i in 0:6]
sols  = Dict()
Xtest = [(x, y) for x in 1.5*λ:0.1:2.5λ, y in -0.75:0.1:-0.25]

for ppw in ppws
    # mesh the obstacles with a smaller mesh size due to their curvature
    h = P * λ / ppw
    meshsize = Dict(ent => h for ent in union(HWW.freesurface(tank), HWW.bottom(tank)))
    foreach(HWW.obstacles(tank)) do ent
        return meshsize[ent] = h / 3
    end
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
ref = 1 ./ ppws[1:end-1] .^ (big(P + 1))
iref = length(ppws) - 1
lines!(ax, ppws[1:end-1], ee[iref] / ref[iref] * ref; label = nothing, linestyle = :dash)
# ylims!(ax,1e-12,1)
xlims!(ax, first(ppws), ppws[end-1])
axislegend(ax; position = :rt)

save("paper/figures/mesh_convergence_jellyfish.pdf", fig)
