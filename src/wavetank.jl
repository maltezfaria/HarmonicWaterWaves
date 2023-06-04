"""
    mutable struct WaveTank
"""
Base.@kwdef mutable struct WaveTank
    # entities
    freesurface::Domain    = Domain()
    bottom::Domain         = Domain()
    obstacles::Domain      = Domain()
    # parameters
    parameters::Parameters   = Parameters()
    # pml change of variables
    pml_func                 = nothing
    # discretized fields
    mesh                     = nothing
    quad                     = nothing
    Sop                      = nothing # single layer operator
    Dop                      = nothing # double layer operator
    S                        = nothing
    D                        = nothing
    f                        = nothing
    σ                        = nothing # |J| ϕ
end

freesurface(tank::WaveTank) = tank.freesurface
bottom(tank::WaveTank)      = tank.bottom
obstacles(tank::WaveTank)   = tank.obstacles
parameters(t::WaveTank)     = t.parameters
pml(t::WaveTank)            = t.pml_func
mesh(t::WaveTank)           = t.mesh
quadrature(t::WaveTank)     = t.quad
ambient_dimension(t::WaveTank) = ambient_dimension(t.freesurface)

function domain(tank::WaveTank)
    freesurface(tank) ∪ bottom(tank) ∪ obstacles(tank)
end

depth(t::WaveTank)        = t |> parameters |> depth
frequency(tank::WaveTank) = tank |> parameters |> frequency
gravity(tank::WaveTank)   = tank |> parameters |> gravity
impedance(tank::WaveTank) = tank |> parameters |> impedance
wavenumber(t::WaveTank)   = t |> parameters |> wavenumber
evanescent_wavenumber(t::WaveTank;kwargs...) = evanescent_wavenumber(parameters(t);kwargs...)
wavelength(t::WaveTank)   = t |> parameters |> wavelength

function set_depth!(wavetank::WaveTank,d)
    p = parameters(wavetank)
    p.depth = d
    return wavetank
end

# helper functions to generate the domain
"""
    add_freesurface!(wavetank,Γf::Domain)
    add_freesurface!(wavetank,lc,hc)

Append `Γf` to the free-surface of `wavetank`. If passed a tuple `(lc,hc)`, then
a plane with low-corner and high-corners given by `lc` and `uc` is created at `z=0`.
"""
function add_freesurface!(wavetank::WaveTank,Γf::Domain)
    union!(freesurface(wavetank),Γf)
    return wavetank
end

function add_freesurface!(wavetank::WaveTank,a::Number,b::Number)
    @assert a < b
    Γf = Domain(line(Point2D(b, 0), Point2D(a, 0)))
    add_freesurface!(wavetank,Γf)
end

function add_freesurface!(wavetank::WaveTank,a::NTuple{2},b::NTuple{2})
    @assert all(a .< b)
    Γf  = ParametricEntity(a,b) do u
        SVector(u[1],u[2],0)
    end |> Domain
    add_freesurface!(wavetank,Γf)
end


"""
    add_bottom!(wavetank,Γb)

Append `Γb` to `bottom(wavetank)`.
"""
function add_bottom!(wavetank, Γb::Domain)
    # TODO: check that endpoints of bottom coincide with depth?
    union!(bottom(wavetank),Γb)
    return wavetank
end
function add_bottom!(wavetank, a::Number,b::Number)
    @assert a < b
    d = depth(wavetank)
    msg = """depth of `wavetank` must be finite.
    Call `set_depth!(tank,d)` to fix the depth"""
    @assert isfinite(d) msg
    Γb = Domain(line(Point2D(a, -d), Point2D(b, -d)))
    add_bottom!(wavetank,Γb)
    return wavetank
end

function add_bottom!(wavetank::WaveTank,a::NTuple{2},b::NTuple{2})
    @assert all(a .< b)
    d = depth(wavetank)
    msg = """depth of `wavetank` must be finite.
    Call `set_depth!(tank,d)` to fix the depth"""
    @assert isfinite(d) msg
    Γf  = ParametricEntity(a,b) do u
        SVector(u[1],u[2],-d)
    end |> Domain
    add_bottom!(wavetank,Γf)
end

function add_pml!(wavetank,pml)
    wavetank.pml_func = pml
end
# function add_pml!(wavetank;a,θ)
#     τ = OrthogonalLinearPML(;a,θ)
#     # τ = OrthogonalQuadraticPML(;a,θ)
#     add_pml!(wavetank,τ)
# end

function add_orthogonal_pml!(wavetank;a,θ)
    τ = OrthogonalLinearPML(;a,θ)
    add_pml!(wavetank,τ)
end

"""
    add_obstacles!(wavetank,Γ::Domain)
"""
function add_obstacles!(wavetank,Γ::Domain)
    union!(obstacles(wavetank),Γ)
    return wavetank
end

function discretize!(tank::WaveTank;meshsize,qorder=9)
    Γ    = domain(tank)
    @info "meshing wavetank"
    mesh = meshgen(Γ;meshsize)
    tank.parameters.meshsize = meshsize
    tank.mesh = mesh
    @info "creating a quadrature for the mesh"
    quad   = NystromMesh(mesh;qorder)
    N = ambient_dimension(mesh)
    tank.parameters.qorder = qorder
    tank.quad = quad
    # assemble lazy integral operators
    @info "creating lazy integral operators"
    p       = parameters(tank)
    τ       = pml(tank)
    if isnothing(τ)
        op      = Laplace(;dim=N)
    else
        op      = LaplacePML(;dim=N,τ)
    end
    Sop     = SingleLayerOperator(op,quad)
    Dop     = DoubleLayerOperator(op,quad)
    tank.Sop = Sop
    tank.Dop = Dop
    return tank
end

function assemble_operators!(tank::WaveTank;correction=:quadgk,compression=:matrix)
    m,n = size(tank.Sop)
    h = tank.parameters.meshsize
    q = tank.parameters.qorder
    @info "assembling $m × $n single- and double-layer operators"
    @info "|--dense computation"
    t = @elapsed begin
        if compression == :matrix
            S       = tank.Sop |> Matrix
            D       = tank.Dop |> Matrix
        end
    end
    @info "    |--- took $t seconds"
    @info "|--sparse correction"
    t = @elapsed begin
        if correction == :quadgk
            δS       = hcubature_correction(tank.Sop;max_dist=5*h,maxevals=2000,atol=min(h^q,1e-12))
            δD       = hcubature_correction(tank.Dop;max_dist=5*h,maxevals=2000,atol=min(h^q,1e-12))
        else
            error("unknown correction method")
        end
    end
    @info "    |--- took $t seconds"
    @info "|--composing dense and sparse operators"
    t = @elapsed begin
        if correction == :none
            tank.S   = S
            tank.D   = D
        elseif compression == :matrix
            tank.S = S + δS
            tank.D = D + δD
        else
            tank.S   = LinearMap(S) + LinearMap(δS)
            tank.D   = LinearMap(D) + LinearMap(δD)
        end
    end
    @info "    |--- took $t seconds"
    return tank
end

function solve!(tank::WaveTank,f::AbstractVector,method=:gmres)
    tank.f  = f
    quad    = quadrature(tank)
    τ       = pml(tank)
    dofs    = quadrature(tank).qnodes
    k       = impedance(tank)
    J_diag  = Diagonal([jacobian_det(τ,dof) for dof in dofs])
    J_diag*f ≈ f || @warn("source term does not vanish inside the PML")
    S,D = tank.S, tank.D
    rhs = S*f
    # compute L = 0.5*|J|^{-1} + D - ω^2/g * Sfs
    L   = 0.5*inv(J_diag) + D
    Γ   = freesurface(tank)
    Is  = dom2qtags(quad,Γ) # index of dofs on free surface
    @. L[:,Is] = L[:,Is] - k*S[:,Is]
    # ϕ   = (L*J_diag)\rhs
    if method === :gmres
        @info "solving with gmres"
        t = @elapsed begin
            σ, hist = gmres(L, rhs; restart=1000, maxiter=1000, log=true, verbose=false, reltol=1e-12)
        end
        @info "|-- took $t seconds"
        @info "|-- gmres converged in $(hist.iters) iterations"
    else
        σ   = L\rhs
    end
    ϕ   = inv(J_diag)*σ
    tank.σ = σ
    return NystromDensity(ϕ,quad)
end

function solve_eigenvalues(tank)
    quad    = quadrature(tank)
    τ       = pml(tank)
    dofs    = quadrature(tank).qnodes
    J_diag  = Diagonal([jacobian_det(τ,dof) for dof in dofs])
    S,D = tank.S, tank.D
    # solve Ax = λBx, with A = 0.5*|J|^{-1} + D
    A   = 0.5*inv(J_diag) + D
    Γ   = freesurface(tank)
    Is  = dom2qtags(quad,Γ) # index of dofs on free surface
    B   = zero(S)
    B[:,Is] = S[:,Is]
    @info "computing generalized eigenvalue decomposition"
    return eigen(A,B)
end

function solve!(tank::WaveTank,f::Function)
    quad = quadrature(tank)
    solve!(tank,[f(dof) for dof in quad.qnodes])
end

function solution(tank)
    G = tank.Sop.kernel
    dG = tank.Dop.kernel
    f = tank.f
    quad    = quadrature(tank)
    τ       = pml(tank)
    dofs    = quadrature(tank).qnodes
    k       = impedance(tank)
    # compute L = 0.5*|J|^{-1} + D - ω^2/g * Sfs
    Γ   = freesurface(tank)
    Is  = dom2qtags(quad,Γ) # index of dofs on free surface
    f′ = deepcopy(f)
    σ  = tank.σ
    @. f′[Is] += k*σ[Is] # ω²/g |J| φ + f
    𝒮 = IntegralPotential(G,quad)[f′]
    𝒟 = IntegralPotential(dG,quad)[σ]
    sol = (x) -> 𝒟(x) - 𝒮(x)
    return sol
end

# incident wave
function plane_wave(tank;θ=0)
    d = depth(tank)
    β = wavenumber(tank)
    A = 1/cosh(β*d) # normalization constant
    return _modal_solution(A,β,d)
end

function _modal_solution(A,β,d)
    ϕᵢ  = (dof) -> begin
        x = coords(dof)
        A*exp(im*β*x[1])*cosh(β*(x[2]+d))
    end
    dϕᵢ = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        A*β*exp(im*β*x[1]) * (im*cosh(β*(x[2]+d))*n[1] + sinh(β*(x[2]+d)*n[2]))
    end
    return ϕᵢ,dϕᵢ
end

function evanescent_wave(tank::WaveTank;n=1,s=1)
    @assert n ≥ 0
    d = depth(tank)
    γ = evanescent_wavenumber(tank;n)
    A = ComplexF64(1) # normalization constant
    return _modal_solution(A,s*im*γ,d)
end

# incident wave
function plane_wave_3d(tank;θ=0)
    d = depth(tank)
    β = wavenumber(tank)
    A = 1/cosh(β*d) # normalization constant
    ϕᵢ  = (dof) -> begin
        x = coords(dof)
        A*exp(im*β*(cos(θ)*x[1] + sin(θ)*x[2]))*cosh(β*(x[3]+d))
    end
    # normal derivative of ϕᵢ in 3d
    dϕᵢ = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        A*β*exp(im*β*(cos(θ)*x[1] + sin(θ)*x[2])) *
        (im*cosh(β*(x[3]+d))*cos(θ)*n[1] + im*cosh(β*(x[3]+d))*sin(θ)*n[2] + sinh(β*(x[3]+d)*n[3]))
    end
    return ϕᵢ,dϕᵢ
end

# plotting

@recipe function f(tank::WaveTank; meshsize=0.1)
    @series begin
        lc := :blue
        label := ""
        view(tank.mesh,tank.freesurface)
    end
    @series begin
        lc := :black
        label := ""
        view(tank.mesh,tank.bottom)
    end
    @series begin
        label := ""
        lc := :red
        meshsize := meshsize
        view(tank.mesh,tank.obstacles)
    end
end
