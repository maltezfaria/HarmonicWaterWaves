using HarmonicWaterWaves
using Test
using StaticArrays
using HCubature
using LinearAlgebra
using Plots
import WavePropBase as WPB

import HarmonicWaterWaves as WW

a, l = 2,4π
h = 0.1
q = 3
θ = π/4

# custom pml on the surface only function
struct SurfacePML
    a::Float64
    c::ComplexF64
end
function SurfacePML(;a,θ)
    SurfacePML(a,exp(im*θ))
end

function (f::SurfacePML)(dof)
    x   = WPB.coords(dof)
    N   = length(x)
    a,c = f.a,f.c
    WPB.svector(N) do d
        xd = x[d]
        if d == N || abs(xd) <= a || !iszero(x[N])
            Complex(xd)
        elseif xd > a
            a + (xd-a)*c
        else
            -a + (xd+a)*c
        end
    end
end

function WW.jacobian_det(f::SurfacePML,dof)
    x = WPB.coords(dof)
    a,c = f.a,f.c
    return (abs(x[1]) > a && x[2] == 0) ? c : one(c)
end

pml = SurfacePML(;a,θ)

p    = HarmonicWaterWaves.Parameters(frequency=√1,gravity=1)
ss = collect(0:0.1:4)
ee = []

for shift in ss
    tank = HarmonicWaterWaves.WaveTank(parameters=p)

    HarmonicWaterWaves.add_freesurface!(tank,-a-l,-a)
    HarmonicWaterWaves.add_freesurface!(tank,-a,a)
    HarmonicWaterWaves.add_freesurface!(tank,a,a+l)

    # HarmonicWaterWaves.add_orthogonal_pml!(tank;a,θ)
    HarmonicWaterWaves.add_pml!(tank,pml)
    x₀  = SVector(shift,-1.0) # center of obstacle
    xₛ  = x₀
    obs = WPB.boundary(WPB.Disk(;radius=0.5,center=x₀,normal=:inside)) |> WPB.Domain
    HarmonicWaterWaves.add_obstacles!(tank,obs)

    HarmonicWaterWaves.discretize!(tank;meshsize=h,qorder=q)

    k = WW.impedance(tank)
    op  = WW.InfiniteDepthWaterWaves(;k)
    G   = WPB.SingleLayerKernel(op)
    dG  = WPB.DoubleLayerKernel(op)
    ϕe   = source -> G(xₛ,source)
    dϕe  = source -> dG(xₛ,source)

    ## Solve with proper Green function
    t_green = @elapsed begin
        Γs         = WW.obstacles(tank)
        quad_obs   = WPB.NystromMesh(view(tank.mesh,Γs);qorder=q)
        op         = WW.InfiniteDepthWaterWaves(;k)
        Sop        = WPB.SingleLayerOperator(op,quad_obs)
        Dop        = WPB.DoubleLayerOperator(op,quad_obs)
        S₀         = Sop |> Matrix
        D₀         = Dop |> Matrix
        δS         = WPB.hcubature_correction(Sop;max_dist=5*h,maxevals=20)
        δD         = WPB.hcubature_correction(Dop;max_dist=5*h,maxevals=20)
        S = S₀ + δS
        D = D₀ + δD
        f          = [dϕe(dof) for dof in quad_obs.qnodes]
        ϕ_green = (I/2 + D)\(S*f)
    end

    xx   = [atan(dof.coords[2]-x₀[2],dof.coords[1]-x₀[1]) for dof in quad_obs.qnodes]
    yy   = [ϕe(dof) for dof in quad_obs.qnodes]

    @info t_green, length(quad_obs.qnodes), norm(yy-ϕ_green)

    scatter(xx,real(yy),label="exact",m=:cross,ms=10)
    scatter!(xx,real(ϕ_green),label="green",m=:circle)
    ##

    ## Solve with PML Green function
    # plot(tank)

    t_pml = @elapsed begin
        quad = tank.quad
        dofs = quad.qnodes
        HarmonicWaterWaves.assemble_operators!(tank;correction=:quadgk)
        f = [dof.coords[2] == 0 ? dϕe(dof) - k*ϕe(dof) : dϕe(dof) for dof in dofs]
        ϕ_pml = HarmonicWaterWaves.solve(tank,f)
    end
    idxs = WPB.dom2qtags(quad,Γs)
    @info t_pml, length(quad.qnodes), norm(yy-ϕ_pml[idxs],Inf)
    push!(ee,norm(yy-ϕ_pml[idxs],Inf))
end

##

# fig = plot(ss,ee,xlabel="center",ylabel="error",
#     label="orthogonal",lw=2,yscale=:log10,m=:circle,legend=:topleft,yticks=[1/(10^i) for i in 5:-1:0])

plot!(fig,ss,ee,label="non-orthogonal",m=:circle,lw=2)
