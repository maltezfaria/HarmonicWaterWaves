using Test
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
HWW.add_obstacles!(tank, obs)

# set a non-uniform mesh size
ppw = 30
P = 5
h = P * λ / ppw
meshsize = Dict(
    ent => h for ent in union(HWW.freesurface(tank), HWW.bottom(tank), HWW.obstacles(tank))
)
foreach(HWW.obstacles(tank)) do ent
    return meshsize[ent] = h / 4
end

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
Ωh = WPB.NystromMesh(tank.obstacles; meshsize = 0.1, qorder = 10) # for identifying points inside the obstacle
step = 0.025
quad = tank.quad
sol = HWW.solution(tank)
xrange = -pml_start-pml_length:step:pml_start+pml_length
yrange = -1.5:step:0
ϕ_sca = NaN * zeros(ComplexF64, length(xrange), length(yrange))
ϕ_tot = NaN * similar(ϕ_sca)
ϕ_inc = NaN * similar(ϕ_sca)
for (i, x) in enumerate(xrange)
    for (j, y) in enumerate(yrange)
        pt = (x, y)
        WPB.isinside(pt, Ωh) && continue
        y < depth_func(x) && continue
        ϕ_sca[i, j] = sol(pt)
        ϕ_inc[i, j] = ϕi(pt)
        ϕ_tot[i, j] = ϕ_sca[i, j] + ϕ_inc[i, j]
    end
end

## plot incident, total, and scattered fields

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
    ymin = 0.33,
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
    ymin = 0.33,
    ymax = 1,
    linestyle = :dash,
    color = (:black, 0.5),
)
lines!(ax3, tank.mesh)

cb = Colorbar(
    fig[1:3, 2];
    limits = (cmin, cmax),
    height = Relative(0.8),
    width = 10,
    # label = L"\Re(\tilde{\varphi})",
)

rowgap!(fig.layout, -30)
colgap!(fig.layout, 10)

save("paper/figures/jellyfish_fields.pdf", fig)

##
if make_animation
    framerate = 30
    timestamps = range(0, 2π; step = 1 / framerate)
    record(
        fig,
        "paper/animations/jellyfish_fields.gif",
        timestamps;
        framerate = framerate,
    ) do t
        return time[] = t
    end
end
