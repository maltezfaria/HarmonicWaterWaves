using HarmonicWaterWaves
using LinearAlgebra
using Plots
using LaTeXStrings
import WavePropBase as WPB
import HarmonicWaterWaves as WW

# gr()

a, l = 1, 1
h    = 0.1
q    = 3
θ    = π/2
d    = 2
b    = a+l

# τ = WW.OrthogonalLinearPML(a,θ,b,d)
τ = WW.OrthogonalLinearPMLStretching(;a,b,θ,d)

default(legendfontsize=12,xlabelfontsize=12,ylabelfontsize=12,guidefontsize=12)

x = -(b+l):0.1:b+l |> collect
y = [τ((x,0))[1] for x in x]
fig = plot(x,real(y),label=L"\mathrm{Re}(\tau)",xlabel=L"x",ls=:solid,lw=2)
plot!(fig,x,imag(y),label=L"\mathrm{Im}(\tau)",lw=2)

xmin, xmax = xlims(fig)
ymin, ymax = ylims(fig)
rectangle(w, h, x, y) = Shape(x .+ [0,w,w,0], y .+ [0,0,h,h])

# plot rectangel in the middle
plot!(fig,rectangle(2a,ymax-ymin,-a,ymin), opacity=.3,label=nothing, color=:blue)

plot!(fig,rectangle(l,ymax-ymin,a,ymin), opacity=.3,label=nothing, color=:green)
plot!(fig,rectangle(l,ymax-ymin,-a-l,ymin), opacity=.3,label=nothing, color=:green)
# vline!(fig,[-a,a,-b,b],ls=:dash,label="",lw=2,lc=:black)

plot!(fig,rectangle(xmax-b,ymax-ymin,b,ymin), opacity=.3,label=nothing, color=:gray)
plot!(fig,rectangle(-b-xmin,ymax-ymin,xmin,ymin), opacity=.3,label=nothing, color=:gray)

xlims!(fig,-3,3)
ylims!(fig,-4,4)
savefig(fig,"paper/figures/pml_real_and_imag.pdf")
fig

# j = [WW.jacobian_det(τ,(x,0)) for x in x]
# fig = plot(x,real(j),label=L"\mathrm{Re}(|\partial \tau|)",xlabel=L"x",ls=:dash,lw=2)
# plot!(fig,x,imag(j),label=L"\mathrm{Im}(|\partial \tau|)",lw=2)

# savefig(fig,"paper/figures/pml_jacobian_det_real_and_imag.pdf")
