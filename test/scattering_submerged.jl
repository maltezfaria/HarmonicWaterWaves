using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 10,5
order = 7             # method order
q     = 2*(order-1)-1 # quadrature order
c = 1
d = 4

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)
HarmonicWaterWaves.set_depth!(tank,d)

λ = WW.wavelength(tank)
h = λ / 10

# create piercing obstacle
r = 1.5# radius of obstacle
x₀  = SVector(0,-2.0) # center of obstacle
obs = let r = r, x₀ = x₀
    WPB.ParametricEntity(0,2π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₀
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

# free surface
HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

# bottom
HarmonicWaterWaves.add_bottom!(tank,-a-l,-a)
HarmonicWaterWaves.add_bottom!(tank,-a,a)
HarmonicWaterWaves.add_bottom!(tank,a,a+l)

# add pml
pml = HarmonicWaterWaves.OrthogonalPML(;a, c)
HarmonicWaterWaves.add_pml!(tank,pml)

# create mesh
HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

# boundary conditions
ϕi, dϕi = HarmonicWaterWaves.plane_wave(tank)

idxs_free   = WPB.dom2qtags(tank.quad,WW.freesurface(tank))
idxs_bottom = WPB.dom2qtags(tank.quad,WW.bottom(tank))
idxs_obs    = WPB.dom2qtags(tank.quad,WW.obstacles(tank))

f = zeros(ComplexF64, length(tank.quad.qnodes))
for i in idxs_obs
    f[i] = -dϕi(tank.quad.qnodes[i])
end

# sanity checks
ee1  = [dϕi(dof) - WW.impedance(tank)*ϕi(dof) for dof in tank.quad.qnodes[idxs_free]]
ee2  = [dϕi(dof) for dof in tank.quad.qnodes[idxs_bottom]]
if norm(ee1,Inf) > 1e-10 || norm(ee2,Inf) > 1e-10
    @warn "Incident wave does not satify surface and/or bottom conditions" norm(ee1,Inf), norm(ee2,Inf)
end

# solve
HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
ϕ_pml = HarmonicWaterWaves.solve!(tank,f)

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
            # uu[i,j] = ϕi(pt) - sol(pt)
            uu[i,j] = -sol(pt)
        else
            uu[i,j] = NaN
        end
    end
end
cmin, cmax = -1,1

anim = @animate for t in 0:0.2:2π
    heatmap(xrange,yrange,real.(exp(-im*t).*uu),aspect_ratio=1,clims=(-1.0,1.0),colorbar=false,size=(800,150))
    plot!(tank,lw=3)
end

gif(anim, "anim_fps15.gif", fps = 15)
