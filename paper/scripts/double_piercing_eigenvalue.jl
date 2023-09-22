using LinearAlgebra
using CairoMakie
import WavePropBase as WPB
import HarmonicWaterWaves as HWW
using Unitful

make_animation = true

ν     = 4
depth = 1 * u"m"
tank  = HWW.WaveTank(; depth, frequency = sqrt(ν) * u"Hz", gravity = 1 * u"m/s^2")

r = 0.25 # radius of obstacles
δ = 2 * r + 1 # distance between center of obstacles

λ = HWW.wavelength(tank)

pml_start  = 2.5
pml_length = 2.5
domain_end = pml_start + pml_length

# add pml
τ = HWW.OrthogonalPML(; pml_start, pml_strength = 1)
HWW.add_pml!(tank, τ)

# piercing obstacles
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

# set a non-uniform mesh size
ppw = 50
P = 5
meshsize = 0.2

# create mesh
HWW.discretize!(tank; meshsize, order = P)
##

idxs_free   = WPB.dom2qtags(tank.quad, HWW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(tank.quad, HWW.bottom(tank))
idxs_obs    = WPB.dom2qtags(tank.quad, HWW.obstacles(tank))

# solve
HWW.assemble_operators!(tank)
F = HWW.solve_eigenvalues(tank)

## plot eigenvalues and an "interesting" eigenfunction

HWW.set_makie_theme!(; width = 469.75499, height = 160)

fig = Figure()
ax1 = Axis(
    fig[1, 1];
    xlabel = L"\mathrm{Re}(\nu)",
    ylabel = L"\mathrm{Im}(\nu)",
    # aspect=DataAspect(),
)
scatter!(ax1, real(F.values), imag(F.values))
xlims!(ax1, 0, 20)
ylims!(ax1, -20, 0)
ax2 = Axis(
    fig[1, 2];
    xlabel = L"\mathrm{Re}(\nu)",
    # aspect=DataAspect(),
)
scatter!(ax2, real(F.values), imag(F.values))
xlims!(ax2, 0, 5)
ylims!(ax2, -5, 0)
I = findall(F.values) do λ
    return 3 < real(λ) < 5 && -3e-1 < imag(λ) < 0
end
i = I[1]
scatter!(
    ax2,
    [real(F.values[i])],
    [imag(F.values[i])];
    markersize = 7,
    color = :red,
    label = "",
)

# now plot the eigenfunction corresponding to the "resonant" eigenvalue
σ = F.vectors[:, I[1]]
tank.σ = σ
tank.f = zero(σ)
step = 0.025
quad = tank.quad
sol = HWW.solution(tank)
xrange = -domain_end:step:domain_end
yrange = -1:step:0
ϕ = NaN * zeros(ComplexF64, length(xrange), length(yrange))
for (i, x) in enumerate(xrange)
    for (j, y) in enumerate(yrange)
        pt = (x, y)
        (norm(pt .- x₁) < r || norm(pt .- x₂) < r) && continue # piercing obstacle
        ϕ[i, j] = sol(pt)
    end
end
ϕmin, ϕmax = extrema(real, filter(!isnan, ϕ))
cmin, cmax = -1, 1

ax1 = Axis(
    fig[2, 1:2];
    xlabel = L"x_1",
    ylabel = L"x_2",
    aspect = DataAspect(),
    xminorticksvisible = false,
    xticksvisible = false,
    yminorticksvisible = false,
    xticklabelsvisible = true,
    topspinevisible = false,
)

h1 = heatmap!(ax1, xrange, yrange, real.(ϕ) / ϕmax; colorrange = (cmin, cmax))

cb = Colorbar(
    fig[2, 3];
    limits = (cmin, cmax),
    height = Relative(1),
    width = 10,
    label = L"\Re(\tilde{\varphi}_\nu)",
)

vlines!(
    ax1,
    [-pml_start, pml_start];
    ymin = 0,
    ymax = 1,
    linestyle = :dash,
    color = (:black, 0.5),
)

rowgap!(fig.layout, 10)
colgap!(fig.layout, 10)

# remesh for plotting the obstacles
HWW.discretize!(tank; meshsize = 0.01, order = :low) # remesh for plotting the obstacles
lines!(ax1, tank.mesh)

save("paper/figures/eigenvalue_problem.pdf", fig)
