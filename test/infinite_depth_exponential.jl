using HarmonicWaterWaves
using LinearAlgebra
using Nystrom
using Nystrom: ParametricSurfaces
using Plots
import WavePropBase as WPB
import HarmonicWaterWaves as WW

plotlyjs()

include(joinpath(HarmonicWaterWaves.PROJECT_ROOT,"test/exact_solution.jl"))

a, l = 10, 4π
h    = 0.1
q    = 3
θ    = π/4

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
tank = HarmonicWaterWaves.WaveTank(parameters=p)

@info "wavelength = $(WW.wavelength(tank))"

HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
HarmonicWaterWaves.add_freesurface!(tank,-a,a)
HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

HarmonicWaterWaves.add_pml!(tank;a,θ)

plot(tank)

HarmonicWaterWaves.discretize!(tank;meshsize=h,order=q)

HarmonicWaterWaves.assemble_operators!(tank)

function f(dof)
    x = dof.coords
    exp(-x[1]^2)
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

idxs = Nystrom.dom2dof(quad,tank.freesurface[2])

#=
    Plot exact solution vs. PML solution vs. proper Green function solution
=#
xx        = [dof.coords[1] for dof in quad.dofs]
@info extrema([dof.coords[1] for dof in quad.dofs[idxs]])

@info norm(ϕ_green[idxs] - ϕ_pml[idxs],Inf) / norm(ϕ_green[idxs],Inf)

fig = plot()
labeled = false
for d in Γ
    idxs = Nystrom.dom2dof(quad,Domain(d))
    plot!(fig,xx[idxs],imag(ϕ_pml.vals[idxs]),label= labeled ? "" : "pml",color=:red,ls=:solid)
    plot!(fig,xx[idxs],imag(ϕ_green[idxs]),label = labeled ? "" : "Green",color=:green,ls=:solid)
    labeled = true
end
fig
