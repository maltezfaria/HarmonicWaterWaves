using HarmonicWaterWaves
using Test
using StaticArrays
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

make_gif = true

atol = 1e-10

d    = 4
ω    = 1
p    = HarmonicWaterWaves.Parameters(frequency=ω,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

λ = WW.wavelength(tank)

# create piercing obstacles
r = 1.5 # radius of obstacle
δ = r + r + 0.5*λ
x₀  = SVector(-δ/2,0.0) # center of obstacle
obs = let r = r, x₀ = x₀
    WPB.ParametricEntity(0,π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₀
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

x₁  = SVector(δ/2,0.0) # center of obstacle
obs = let r = r, x₁ = x₁
    WPB.ParametricEntity(0,π) do (u,)
        SVector(r*cos(u),-r*sin(u)) + x₁
    end
end |> WPB.Domain
HarmonicWaterWaves.add_obstacles!(tank,obs)

L⁻ = x₀[1] - r
L⁺ = x₁[1] + r

a = max(abs(L⁻),abs(L⁺)) + λ/10
l = 2*λ
c = 1
h = λ/10
q = 11

# free surface
HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,x₀[1]-r)
HarmonicWaterWaves.add_freesurface!(tank,x₀[1]+r,x₁[1]-r)
HarmonicWaterWaves.add_freesurface!(tank,x₁[1]+r,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

# bottom
HarmonicWaterWaves.set_depth!(tank,d)
HarmonicWaterWaves.add_bottom!(tank,-a-l,-a)
HarmonicWaterWaves.add_bottom!(tank,-a,a)
HarmonicWaterWaves.add_bottom!(tank,a,a+l)

# add pml
pml = HarmonicWaterWaves.OrthogonalPML(;a, c)
HarmonicWaterWaves.add_pml!(tank,pml)

# create mesh
HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

# solve
HarmonicWaterWaves.assemble_operators!(tank;atol)
F = HarmonicWaterWaves.solve_eigenvalues(tank)

## plot eigenvalues and an "interesting" eigenfunction
fig = scatter(F.values,label="")
xlims!(fig,0,5)
ylims!(fig,-5,0)
I = findall(F.values) do λ
    1 < real(λ) < 5 && -1e-1 < imag(λ) < 0
end
i = I[1]
scatter!(fig,[real(F.values[i])],[imag(F.values[i])],ms=7,mc=:red,m=:x,label="")
savefig(fig,"paper/figures/eigenvalues_two_piercing_depth_$d.pdf")

## plot eigenfunction on the surface
σ = F.vectors[:,I[1]]
tank.σ = σ
tank.f = zero(σ)

Ω₀ =  WPB.Disk(;center=x₀,radius=1.0*r) # for exluding points inside the obstacle
Ω₁ =  WPB.Disk(;center=x₁,radius=1.0*r) # for exluding points inside the obstacle
step = 0.05
quad = tank.quad
sol = WW.solution(tank)
xrange = -a-l:step:a+l
yrange = -d:step:0
ϕ_eig = NaN*zeros(ComplexF64,length(yrange),length(xrange))
for (i,y) in enumerate(yrange)
    for (j,x) in enumerate(xrange)
        pt = SVector(x,y)
       ( pt ∈ Ω₀ || pt ∈ Ω₁) && continue
        ϕ_eig[i,j] = sol(pt)
    end
end
cmin, cmax = -0.2,0.2
fig_eig = heatmap(xrange,yrange,real.(ϕ_eig),aspect_ratio=1,colorbar=true,size=(800,150),xlims=(-a-l,a+l),ylims=(-d,0),clims=(cmin,cmax))
plot!(fig_eig,tank,lw=3)
vline!(fig_eig,[-a,a],label=nothing)
savefig(fig_eig,"paper/figures/eigenfunction_two_piercing_depth_$d.pdf")
fig_eig
