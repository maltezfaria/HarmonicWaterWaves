using LinearAlgebra
using CairoMakie
import WavePropBase as WPB
import HarmonicWaterWaves as HWW
using Unitful

make_animation = true

ν     = 4
depth = 1 * u"m"
tank  = HWW.WaveTank(; depth, frequency = sqrt(ν) * u"Hz", gravity = 1 * u"m/s^2")

λ = HWW.wavelength(tank)

pml_start  = 3λ
pml_length = λ
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

# set a non-uniform mesh size
ppw = 30
P = 5
meshsize = P * λ / ppw

# create mesh
HWW.discretize!(tank; meshsize, order = P)
##

# boundary conditions
ϕi, dϕi = HWW.plane_wave(tank)

idxs_free   = WPB.dom2qtags(tank.quad, HWW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(tank.quad, HWW.bottom(tank))
idxs_obs    = WPB.dom2qtags(tank.quad, HWW.obstacles(tank))

f = zeros(ComplexF64, length(tank.quad.qnodes))
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

# compute solution on volume
step = 0.025
quad = tank.quad
sol = HWW.solution(tank)
xrange = -pml_start-pml_length:step:pml_start+pml_length
yrange = -1:step:0
ϕ_sca = NaN * zeros(ComplexF64, length(xrange), length(yrange))
ϕ_tot = NaN * similar(ϕ_sca)
ϕ_inc = NaN * similar(ϕ_sca)
for (i, x) in enumerate(xrange)
    for (j, y) in enumerate(yrange)
        pt = (x, y)
        (norm(pt .- x₁) < r || norm(pt .- x₂) < r) && continue # piercing obstacle
        ϕ_inc[i, j] = ϕi(pt)
        ϕ_sca[i, j] = sol(pt)
        ϕ_tot[i, j] = ϕ_sca[i, j] + ϕ_inc[i, j]
    end
end
cmin, cmax = -1, 1

## plot the scattered field

HWW.discretize!(tank; meshsize = 0.01, order = :low) # remesh for plotting the obstacles
time = Observable(0.0) # all plots will be updated when this changes

cmin, cmax = -1, 1

HWW.set_makie_theme!(; width = 469.75499, height = 240)

fig = Figure()

ax1 = Axis(
    fig[1, 1];
    ylabel = L"x_2",
    aspect = DataAspect(),
    xminorticksvisible = true,
    xticksvisible = true,
    yminorticksvisible = false,
    xticklabelsvisible = false,
    topspinevisible = false,
)

h1 = heatmap!(
    ax1,
    xrange,
    yrange,
    lift(t -> real.(exp(-im * t) * ϕ_inc), time);
    colorrange = (cmin, cmax),
)
lines!(ax1, tank.mesh)

ax2 = Axis(
    fig[2, 1];
    ylabel = L"x_2",
    aspect = DataAspect(),
    yminorticksvisible = false,
    topspinevisible = false,
    xticklabelsvisible = false,
)
h2 = heatmap!(
    ax2,
    xrange,
    yrange,
    lift(t -> real.(exp(-im * t) * ϕ_tot), time);
    colorrange = (cmin, cmax),
)
vlines!(
    ax2,
    [-pml_start, pml_start];
    ymin = 0,
    ymax = 1,
    linestyle = :dash,
    color = (:black, 0.5),
)
lines!(ax2, tank.mesh)

ax3 = Axis(
    fig[3, 1];
    xlabel = L"x_1",
    ylabel = L"x_2",
    aspect = DataAspect(),
    xminorticksvisible = true,
    xticksvisible = true,
    yminorticksvisible = false,
    xticklabelsvisible = true,
    topspinevisible = false,
)

h3 = heatmap!(
    ax3,
    xrange,
    yrange,
    lift(t -> real.(exp(-im * t) * ϕ_sca), time);
    colorrange = (cmin, cmax),
)
vlines!(
    ax3,
    [-pml_start, pml_start];
    ymin = 0,
    ymax = 1,
    linestyle = :dash,
    color = (:black, 0.5),
)
lines!(ax3, tank.mesh)

cb = Colorbar(
    fig[1:3, 2];
    limits = (cmin, cmax),
    height = Relative(0.6),
    width = 10,
    label = L"\Re(\tilde{\varphi})",
)

rowgap!(fig.layout, -70)
colgap!(fig.layout, 10)

save("paper/figures/double_piercing_fields.pdf", fig)

##
if make_animation
    framerate = 30
    timestamps = range(0, 2π; step = 1 / framerate)
    record(
        fig,
        "paper/animations/double_piercing_animation.gif",
        timestamps;
        framerate = framerate,
    ) do t
        return time[] = t
    end
end
