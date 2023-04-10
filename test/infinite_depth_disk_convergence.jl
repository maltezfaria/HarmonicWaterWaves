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

ee = []
ll = [π,2π,4π,8π]
for l in ll
    a = 4
    h = 0.1
    q = 5
    θ = π/3

    p    = HarmonicWaterWaves.Parameters(frequency=sqrt(1),gravity=1)
    tank = HarmonicWaterWaves.WaveTank(parameters=p)

    HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
    HarmonicWaterWaves.add_freesurface!(tank,-a,a)
    HarmonicWaterWaves.add_freesurface!(tank,a,a+l)
    # HarmonicWaterWaves.add_freesurface!(tank,-a-l,a+l)

    HarmonicWaterWaves.add_pml!(tank;a,θ)

    x₀  = SVector(1,-1.0) # center of obstacle
    # xₛ  = SVector(0,-1.0) # center of source
    xₛ  = x₀
    obs = boundary(ParametricSurfaces.Disk(;radius=0.5,center=x₀,normal=:inside))
    HarmonicWaterWaves.add_obstacles!(tank,obs)

    HarmonicWaterWaves.discretize!(tank;meshsize=h,order=q)

    k = WW.impedance(tank)
    @info WW.wavelength(tank)
    op  = WW.InfiniteDepthWaterWaves(;k)
    G   = SingleLayerKernel(op)
    dG  = DoubleLayerKernel(op)
    ϕe   = source -> G(xₛ,source)
    dϕe  = source -> dG(xₛ,source)

    ## Solve with PML Green function
    # plot(tank)
    quad = tank.quad
    dofs = quad.dofs
    HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
    f = [dof.coords[2] == 0 ? dϕe(dof) - k*ϕe(dof) : dϕe(dof) for dof in dofs]
    ϕ_pml = HarmonicWaterWaves.solve(tank,f)

    Γs   = WW.obstacles(tank)
    idxs = Nystrom.dom2dof(quad,Γs)

    yy   = [ϕe(dof) for dof in quad.dofs]
    push!(ee,norm(yy[idxs]-ϕ_pml[idxs],Inf))
end

plot(ll,ee,yscale=:log10,m=:cross)
