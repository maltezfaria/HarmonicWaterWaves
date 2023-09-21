using LinearAlgebra
using CairoMakie
using LaTeXStrings
using Unitful
import WavePropBase as WPB
import HarmonicWaterWaves as HWW

HWW.set_makie_theme!()

depth = 1u"m"
ν     = 2
tank  = HWW.WaveTank(; frequency = sqrt(ν)*u"Hz", depth, gravity = 1*u"m/s^2")
λ = HWW.wavelength(tank)
pml_start  = 2*λ
pml_length = 4λ
domain_end = pml_start + pml_length
τ = HWW.OrthogonalPML(;pml_start)

HWW.add_pml!(tank,τ)
HWW.add_wavemaker!(tank,domain_end)

ppws = [4*2^i for i in 0:6]
PP   = [2 5 10]
dict = Dict(p=>Float64[] for p in PP)

for P in PP
    for ppw in ppws
        h = P*λ / ppw
        HWW.discretize!(tank;meshsize=h,order=P)
        HWW.assemble_operators!(tank)
        ϕ, dϕ = HWW.plane_wave(tank)
        ψ, dψ = HWW.evanescent_wave(tank;s=1)
        ϕi = (dof)  -> ϕ(dof) + ψ(dof)
        dϕi = (dof) -> dϕ(dof) + dψ(dof)
        quad = tank.quad
        idxs_free   = WPB.dom2qtags(quad,HWW.freesurface(tank))
        idxs_bottom = WPB.dom2qtags(quad,HWW.bottom(tank))
        idxs_obs   = WPB.dom2qtags(quad,HWW.obstacles(tank))
        f = zeros(ComplexF64, length(tank.quad.qnodes))
        for i in idxs_obs
            f[i] = dϕi(tank.quad.qnodes[i])
        end
        ϕ_pml = HWW.solve!(tank,f)
        ϕ_exact = [ϕi((dof)) for dof in quad.qnodes]
        x̃ = [τ(dof) for dof in quad.qnodes]
        ϕ_pml_exact = [ϕi(x) for x in x̃]
        # error in free surface
        push!(dict[P],norm(ϕ_pml_exact - ϕ_pml,Inf))
    end
end

##
fig = Figure()
xlabels = [string(Int(p)) for p in ppws]
ax  = Axis(fig[1,1],
    xscale=log10,
    yscale=log10,
    xticks=ppws,
    xminorticksvisible=false,
    yminorticksvisible=false,
    xlabel="points per wavelength",
    ylabel="error"
)

for (i,P) in enumerate(PP)
    vals = dict[P]
    scatterlines!(ax,ppws,vals,label="P = $P")
    ref   = 1 ./ ppws.^(big(P))
    iref = P == 10 ? 4 : length(ppws) - 1
    lines!(ax,ppws,vals[iref]/ref[iref]*ref,label=nothing, linestyle=:dash)
end
ylims!(ax,1e-12,1)
axislegend(ax,position=:lb)

save("paper/figures/mesh_convergence_waveguide.pdf",fig)
