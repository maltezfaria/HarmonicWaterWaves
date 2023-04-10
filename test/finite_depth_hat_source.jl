using HarmonicWaterWaves
using LinearAlgebra
using Nystrom
using Nystrom: ParametricSurfaces
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW

plotlyjs()

include(joinpath(HarmonicWaterWaves.PROJECT_ROOT,"test/exact_solution.jl"))

a, l = 10, 8π
h    = 0.1
q    = 3
θ    = π/4
d    = 1

SIDES  = false
BOTTOM = true

p    = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

BOTTOM && HarmonicWaterWaves.add_bottom!(tank,-a-l,a+l,d)

if SIDES
    s1 = ParametricSurfaces.line(WW.Point2D(a+l,-d),WW.Point2D(a+l,0))
    s2 = ParametricSurfaces.line(WW.Point2D(-a-l,0),WW.Point2D(-a-l,-d))
    HarmonicWaterWaves.add_obstacles!(tank,Domain([s1,s2]))
end

HarmonicWaterWaves.add_pml!(tank;a,θ)

plot(tank)

HarmonicWaterWaves.discretize!(tank;meshsize=h,order=q)

HarmonicWaterWaves.assemble_operators!(tank)

function f(dof)
    x = dof.coords
    if -1 < x[1] < 0
        return x[1] + 1
    elseif 0 ≤ x[1] < 1
        return -x[1] + 1
    else
        return zero(x[1])
    end
end

ϕ_pml = HarmonicWaterWaves.solve(tank,f)

#=
    Plot exact solution vs. PML solution vs. proper Green function solution
=#
quad = tank.quad
k    = WW.impedance(tank)
xx        = [dof.coords[1] for dof in quad.dofs]
@info extrema([dof.coords[1] for dof in quad.dofs[idxs]])
ϕ_exact   = solution_totale.(xx,0.0,k)

@info norm(real(ϕ_exact[idxs] - ϕ_pml[idxs]),Inf)

fig = plot()
labeled = false
Γ    = WW.freesurface(tank)
for d in Γ
    idxs = Nystrom.dom2dof(quad,Domain(d))
    plot!(fig,xx[idxs],real(ϕ_pml.vals[idxs]),label= labeled ? "" : "pml",color=:red,ls=:solid)
    plot!(fig,xx[idxs],ϕ_exact[idxs],label= labeled ? "" : "exact",color=:blue,lw=2,ls=:dash)
    labeled = true
end
fig
