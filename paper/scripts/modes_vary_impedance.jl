import HarmonicWaterWaves as WW
using Unitful
using CairoMakie

WW.set_makie_theme!()

frequency    = sqrt(9.81)u"Hz"

kk = Float64[]
γγ = Float64[]
νν = 0.001:0.01:8

for ν in νν
    depth = ν*u"m"
    p    = WW.WaveTank(;frequency,depth,gravity=9.81u"m/s^2")
    @assert WW.impedance(p) ≈ ν
    k    = WW.wavenumber(p)
    γ    = WW.evanescent_wavenumber(p)
    push!(kk,k)
    push!(γγ,γ)
end

##

fig = Figure()
ax  = Axis(fig[1,1],xlabel=L"\nu")
lines!(ax,νν,kk,label=L"k")
lines!(ax,νν,γγ,label=L"\gamma_1")
xlims!(0,8)
# ylims!(0,4)
axislegend(ax; position=:lt)

save("paper/figures/propagative_vs_evanescent.pdf",fig)

fig
