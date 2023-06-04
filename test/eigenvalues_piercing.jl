using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 10,5
h = 0.1
q = 3
θ = π/4
d = 2

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

# create piercing obstacle
r = 1.5 # radius of obstacle
x₀  = SVector(0,0.0) # center of obstacle
obs = let r = r, x₀ = x₀
    WPB.ParametricEntity(0,π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₀
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

# free surface
HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,-r)
# HarmonicWaterWaves.add_freesurface!(tank,-r,r)
HarmonicWaterWaves.add_freesurface!(tank,r,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

# bottom
HarmonicWaterWaves.set_depth!(tank,d)
HarmonicWaterWaves.add_bottom!(tank,-a-l,-a)
HarmonicWaterWaves.add_bottom!(tank,-a,a)
HarmonicWaterWaves.add_bottom!(tank,a,a+l)

# add pml
HarmonicWaterWaves.add_orthogonal_pml!(tank;a,θ)

# create mesh
HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

idxs_free   = WPB.dom2qtags(tank.quad,WW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(tank.quad,WW.bottom(tank))
idxs_obs   = WPB.dom2qtags(tank.quad,WW.obstacles(tank))

# solve
HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
F = HarmonicWaterWaves.solve_eigenvalues(tank)

##
scale_size = 1.0
Plots.scalefontsizes(scale_size)
Γf = WW.freesurface(tank)
xx = [dof.coords[1] for dof in tank.quad.qnodes]
# fig = plot(tank,aspect_ratio=4)
fig = plot(;legendfontsize=12)
labeled = false
for Γ in Γf
    I = WPB.dom2qtags(tank.quad,Γ)
    plot!(fig,xx[I],real(ϕi.(tank.quad.qnodes[I]) +  ϕ_pml.vals[I]),label= labeled ? "" : "approximate pml",color=:red,ls=:solid)
    labeled = true
end
plot!(fig,tank)

##
quad = tank.quad
sol = WW.solution(tank)
xrange = -a-l:0.1:a+l
yrange = -d:0.1:0
uu = zeros(ComplexF64,length(yrange),length(xrange))
for (i,y) in enumerate(yrange)
    for (j,x) in enumerate(xrange)
        pt = SVector(x,y)
        if WPB.isinside(pt,tank.quad)
            uu[i,j] = ϕi(pt) - sol(pt)
        else
            uu[i,j] = NaN
        end
    end
end
cmin, cmax = -1,1

anim = @animate for t in 0:0.2:2π
    heatmap(xrange,yrange,real.(exp(-im*t).*uu),aspect_ratio=1,clims=(cmin,cmax),colorbar=false,size=(800,150),xlims=(-a,a))
    plot!(tank,lw=3)
end

gif(anim, "anim_total_r_$r.gif", fps = 15)
