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
end

freesurface(tank::WaveTank) = tank.freesurface
bottom(tank::WaveTank)      = tank.bottom
obstacles(tank::WaveTank)   = tank.obstacles
parameters(t::WaveTank)     = t.parameters
pml(t::WaveTank)            = t.pml_func
mesh(t::WaveTank)           = t.mesh
quadrature(t::WaveTank)     = t.quad

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

Append `Γf` to the free-surface of `wavetank`.
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

function discretize!(tank::WaveTank;meshsize,qorder)
    Γ    = domain(tank)
    @info "meshing wavetank"
    mesh = meshgen(Γ;meshsize)
    tank.parameters.meshsize = meshsize
    tank.mesh = mesh
    @info "creating a quadrature for the mesh"
    quad   = NystromMesh(mesh;qorder)
    tank.parameters.qorder = qorder
    tank.quad = quad
    # assemble lazy integral operators
    @info "creating lazy integral operators"
    p       = parameters(tank)
    τ       = pml(tank)
    if isnothing(τ)
        op      = Laplace(;dim=2)
    else
        op      = LaplacePML(;dim=2,τ)
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
    @info "assembling $m × $n single- and double-layer operators"
    if compression == :matrix
        S       = tank.Sop |> Matrix
        D       = tank.Dop |> Matrix
    end
    if correction == :quadgk
        δS       = hcubature_correction(tank.Sop;max_dist=10*h,maxevals=200)
        δD       = hcubature_correction(tank.Dop;max_dist=10*h,maxevals=200)
    end
    tank.S = S + δS
    tank.D = D + δD
    return tank
end

function solve(tank::WaveTank,f::AbstractVector)
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
    σ   = L\rhs
    ϕ   = inv(J_diag)*σ
    return NystromDensity(ϕ,quad)
end

function solve(tank::WaveTank,f::Function)
    quad = quadrature(tank)
    solve(tank,[f(dof) for dof in quad.qnodes])
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
    return _modal_solution(A,im*γ,d)
    # ϕᵢ  = (dof) -> begin
    #     x = coords(dof)
    #     A*exp(-γ*x[1])*cos(γ*(x[2]+d))
    # end
    # dϕᵢ = (dof) -> begin
    #     x = coords(dof)
    #     n = normal(dof)
    #     -A*γ*exp(-γ*x[1]) * (cos(γ*(x[2]+d))*n[1] - sin(γ*(x[2]+d)*n[2]))
    # end
    # return ϕᵢ,dϕᵢ
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
