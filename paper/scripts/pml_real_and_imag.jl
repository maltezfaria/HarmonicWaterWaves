using HarmonicWaterWaves
using LinearAlgebra
using CairoMakie
using LaTeXStrings
import WavePropBase as WPB
import HarmonicWaterWaves as WW

th = WW.set_makie_theme!()

a, b, c, d = 1, 2, 1, 2
τ = WW.OrthogonalPML(;pml_start=a,stretch_start=b,pml_strength=c,stretch_strength=d)

x = -(b+1):0.1:b+1 |> collect
y = [τ((x,0))[1] for x in x]

fig, ax, l1 = lines(x,real(y);
    label=L"\mathrm{Re}(\tau)",
    axis = (xlabel=L"x",
            limits=(-3,3,-4,4),
            # xticks=-3:1:3,
            # yticks=-4:1:4
            ),
)
l2 = lines!(ax,x,imag(y);
    label=L"\mathrm{Im}(\tau)"
)
# fix limits
# autolimits!(ax)
# xlims!(ax,-3,3)
# ylims!(ax,-4,4)

# shade the three regions
vspan!(ax,-1, 1, color = (:gray, 0.2))
vspan!(ax,-2, -1, color = (:blue, 0.2))
vspan!(ax,1, 2, color = (:blue, 0.2))
vspan!(ax,-3, -2, color = (:red, 0.2))
vspan!(ax,2, 3, color = (:red, 0.2))

# add a legend
leg = axislegend(ax;
    position=:lt,
    framevisible=true,
    padding=(1,1,-4,-4)
)

# resize_to_layout!(fig)

save("paper/figures/pml_real_and_imag.pdf",fig)
