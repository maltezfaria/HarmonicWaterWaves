using HarmonicWaterWaves
using Test
using StaticArrays
using HCubature
using LinearAlgebra
using Nystrom
using Nystrom: ParametricSurfaces
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW

plotlyjs()

include(joinpath(HarmonicWaterWaves.PROJECT_ROOT,"test/exact_solution.jl"))

a, l = 10, 10
h = 0.1
q = 3

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1,θ=π/4)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

HarmonicWaterWaves.add_free_surface!(tank,a)
HarmonicWaterWaves.add_pml_layer!(tank,l)

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
    Solve it with proper Green function
=#

k  = HarmonicWaterWaves.impedance(tank)
Γ = HarmonicWaterWaves.domain(tank)
quad   = tank.quad
op = WW.InfiniteDepthWaterWaves(;k)
S  = SingleLayerOperator(op,quad) |> Nystrom.assemble_gk

# note: there is no 1/2 here because of the behaviour of the Green function on
# the freesurface; see equation 8.2 of Hein, Duran, Nedelec
ϕ_green = S*[f(dof) for dof in quad.dofs]



#=
    Plot exact solution vs. PML solution vs. proper Green function solution
=#
xx        = [dof.coords[1] for dof in quad.dofs]
ϕ_exact   = solution_totale.(xx,0.0,k)

fig = plot()
labeled = false
for d in Γ
    idxs = Nystrom.dom2dof(quad,Domain(d))
    plot!(fig,xx[idxs],real(ϕ_pml.vals[idxs]),label= labeled ? "" : "pml",color=:red,ls=:solid)
    plot!(fig,xx[idxs],real(ϕ_green[idxs]),label = labeled ? "" : "Green",color=:green,ls=:solid)
    plot!(fig,xx[idxs],ϕ_exact[idxs],label= labeled ? "" : "exact",color=:blue,lw=2,ls=:dash)
    labeled = true
end
fig
