"""
    mutable struct WaveTank

Main type exported by `HarmonicWaterWaves` describing a numerical experiment.
"""
Base.@kwdef mutable struct WaveTank
    # entities
    freesurface::Domain = Domain()
    bottom::Domain      = Domain()
    obstacles::Domain   = Domain()
    # Dimensional parameters
    frequency::typeof(1.0u"Hz")
    gravity::typeof(1.0u"m/s^2") = 9.81 * u"m/s^2"
    depth::typeof(1.0u"m")
    lengthscale::typeof(1.0u"m") = depth
    # mesh parameters
    meshsize::Union{Float64,Dict} = Inf
    qorder::Int       = 0
    # pml change of variables
    pml_func = nothing
    # discretized fields
    mesh = nothing
    quad = nothing
    Sop  = nothing # single layer operator
    Dop  = nothing # double layer operator
    S    = nothing
    D    = nothing
    f    = nothing
    σ    = nothing # |J| ϕ
end

frequency(p::WaveTank) = p.frequency
gravity(p::WaveTank) = p.gravity
depth(p::WaveTank) = p.depth
lengthscale(p::WaveTank) = p.lengthscale

function Base.empty!(tank::WaveTank)
    (tank.freesurface = Domain(); tank.bottom = Domain(); tank.obstacles = Domain())
end

"""
    impendance(p::WaveTank)

Impedance-like dimensioness parameter given by `x₀*ω^2/g`, where `x₀` is the
characteristic (dimensional) spatial scale, `ω` is the frequency, and `g` is the
gravity.
"""
function impedance(p::WaveTank)
    g = gravity(p)
    λ = lengthscale(p)
    ω = frequency(p)
    return ω^2 * λ / g |> NoUnits
end

"""
    wavenumber(p::WaveTank)

Dimensionless wavenumber of the propagating mode in the wave tank `p`.

To obtain the dimensional wavenumber, divide the result by `lengthscale(p)`.
"""
function wavenumber(p::WaveTank)
    ν = impedance(p)
    if isinf(depth(p))
        error("Infinite depth not (yet) supported")
    else
        k = find_zero(k -> k * tanh(k) - ν, ν)
    end
    return k
end

"""
    wavelength(p::WaveTank)

Dimensionless wavelength of the propagating mode in the wave tank
`p`.

To obtain the dimensional wavelength, multiply the result by `lengthscale(p)`.
"""
wavelength(p::WaveTank) = 2π / wavenumber(p)

set_lengthscale!(tank::WaveTank, λ = wavelength(tank)) = (tank.lengthscale = λ)

"""
    evanescent_wavenumber(p::WaveTank; n = 1)

The `n`-th evanescent (dimensionless) wavenumber of the wave tank `p`.
"""
function evanescent_wavenumber(p::WaveTank; n = 1)
    @assert n ≥ 1
    ν = impedance(p)
    if isinf(depth(p))
        error("Infinite depth not (yet) supported")
    else
        f = γ -> γ * tan(γ) + ν
        a = (1 + eps()) * (n * π - π / 2)
        b = (1 - eps()) * (n * π)
        @debug f(a), f(b)
        γ = find_zero(f, (a, b))
    end
    return γ
end

"""
    evanescent_wavelength(p::WaveTank; n = 1)

The `n`-th evanescent (dimensionless) "wavelength" of the wave tank `p`.
"""
evanescent_wavelength(p::WaveTank; n = 1) = 2π / evanescent_wavenumber(p; n)

freesurface(tank::WaveTank)    = tank.freesurface
bottom(tank::WaveTank)         = tank.bottom
obstacles(tank::WaveTank)      = tank.obstacles
pml(t::WaveTank)               = t.pml_func
mesh(t::WaveTank)              = t.mesh
quadrature(t::WaveTank)        = t.quad
ambient_dimension(t::WaveTank) = ambient_dimension(t.freesurface)

function domain(tank::WaveTank)
    return freesurface(tank) ∪ bottom(tank) ∪ obstacles(tank)
end

# helper functions to generate the domain
"""
    add_freesurface!(wavetank,Γf::Domain)
    add_freesurface!(wavetank,lc,hc)

Append `Γf` to the free-surface of `wavetank`. If passed a tuple `(lc,hc)`, then
a plane with low-corner and high-corners given by `lc` and `uc` is created at
`z=0`.

The normal vector should be oriented downwards into the fluid domain.
"""
function add_freesurface!(wavetank::WaveTank, Γf::Domain)
    union!(freesurface(wavetank), Γf)
    return wavetank
end

function add_freesurface!(wavetank::WaveTank, a::Number, b::Number)
    @assert a < b
    Γf = Domain(line(Point2D(a, 0), Point2D(b, 0)))
    return add_freesurface!(wavetank, Γf)
end

function add_freesurface!(wavetank::WaveTank, a::NTuple{2}, b::NTuple{2})
    @assert all(a .< b)
    Γf = ParametricEntity(a, b) do u
        return SVector(u[1], u[2], 0)
    end |> Domain
    return add_freesurface!(wavetank, Γf)
end

"""
    add_bottom!(wavetank,Γb)
    add_bottom!(wavetank,a::Number,b::Number; depth)

Append `Γb` to `bottom(wavetank)`. If passed two numbers `a` and `b`, create a
line connecting `(a,-d)` to `(b,-d)`, where `d` is the depth of `wavetank`.
"""
function add_bottom!(wavetank, Γb::Domain)
    # TODO: check that endpoints of bottom coincide with depth?
    union!(bottom(wavetank), Γb)
    return wavetank
end

function add_bottom!(wavetank, a::Number, b::Number; depth = 1)
    @assert a < b
    msg = """depth of `wavetank` must be finite."""
    @assert isfinite(depth) msg
    Γb = Domain(line(Point2D(b, -depth), Point2D(a, -depth)))
    add_bottom!(wavetank, Γb)
    return wavetank
end

"""
    add_pml!(wavetank,pml)

Add a PML to the wavetank. The PML is a function that takes a point `x ∈ ℝᵈ` and
returns a point `x̃ ∈ ℂᵈ`.

See also [`OrthogonalPML`](@ref).
"""
function add_pml!(wavetank, pml)
    return wavetank.pml_func = pml
end

"""
    add_obstacles!(wavetank,Γ::Domain)

Add an obstacle to `wavetank`.
"""
function add_obstacles!(wavetank, Γ::Domain)
    union!(obstacles(wavetank), Γ)
    return wavetank
end

function discretize!(tank::WaveTank; meshsize, order = :medium)
    if order === :low
        P = 2
    elseif order === :medium
        P = 5
    elseif order === :high
        P = 10
    elseif isa(order, Number)
        P = order
    else
        error("order must be either an integer, or one of :low, :medium, :high")
    end
    qorder = 2P - 1
    Γ = domain(tank)
    @info "meshing wavetank"
    mesh = meshgen(Γ; meshsize)
    tank.meshsize = meshsize
    tank.mesh = mesh
    @info "creating a quadrature for the mesh"
    quad = NystromMesh(mesh; qorder)
    N = ambient_dimension(mesh)
    tank.qorder = qorder
    tank.quad = quad
    # assemble lazy integral operators
    @info "creating lazy integral operators"
    τ = pml(tank)
    op = LaplacePML(; dim = N, τ)
    Sop = SingleLayerOperator(op, quad)
    Dop = DoubleLayerOperator(op, quad)
    tank.Sop = Sop
    tank.Dop = Dop
    return tank
end

function assemble_operators!(
    tank::WaveTank;
    compression = :matrix,
    atol = 1e-10,
    maxevals = 2000,
    initdiv = 1,
    maxdist = nothing,
)
    m, n = size(tank.Sop)
    h = maximum(values(tank.meshsize))
    isnothing(maxdist) && (maxdist = 5 * h)
    @info "assembling $m × $n single- and double-layer operators"
    @info "|--dense computation"
    t = @elapsed begin
        if compression == :matrix
            S = tank.Sop |> Matrix
            D = tank.Dop |> Matrix
        end
    end
    @info "    |--- took $t seconds"
    @info "|--sparse correction"
    t = @elapsed begin
        δS = hcubature_correction(tank.Sop; maxdist, atol, maxevals, initdiv)
        δD = hcubature_correction(tank.Dop; maxdist, atol, maxevals, initdiv)
    end
    @info "    |--- took $t seconds"
    @info "|--composing dense and sparse operators"
    t = @elapsed begin
        if compression == :matrix
            tank.S = S + δS
            tank.D = D + δD
            # tank.S = S
            # tank.D = D
        else
            tank.S = LinearMap(S) + LinearMap(δS)
            tank.D = LinearMap(D) + LinearMap(δD)
        end
    end
    @info "    |--- took $t seconds"
    return tank
end

function solve!(tank::WaveTank, f::AbstractVector, method = :gmres)
    tank.f = f
    quad   = quadrature(tank)
    τ      = pml(tank)
    dofs   = quadrature(tank).qnodes
    ν      = impedance(tank)
    J_diag = Diagonal([jacobian_det(τ, dof) for dof ∈ dofs])
    J_diag * f ≈ f || @warn("source term does not vanish inside the PML")
    S, D = tank.S, tank.D
    rhs = S * (J_diag * f)
    # compute L = -0.5*|J|^{-1} + D + ω^2/g * Sfs
    L = -0.5 * inv(J_diag) + D
    Γ = freesurface(tank)
    Is = dom2qtags(quad, Γ) # index of dofs on free surface
    @. L[:, Is] = L[:, Is] + ν * S[:, Is]
    # ϕ   = (L*J_diag)\rhs
    if method === :gmres
        @info "solving with gmres"
        t = @elapsed begin
            σ, hist = gmres(
                L,
                rhs;
                restart = 1000,
                maxiter = 1000,
                log = true,
                verbose = false,
                reltol = 1e-12,
            )
        end
        @info "|-- took $t seconds"
        @info "|-- gmres converged in $(hist.iters) iterations"
    elseif method === :direct
        @info "solving with direct solver"
        t = @elapsed begin
            σ = L \ rhs
        end
        @info "|-- took $t seconds"
    else
        error("unknown method")
    end
    ϕ = inv(J_diag) * σ
    tank.σ = σ
    return NystromDensity(ϕ, quad)
end

function solve_eigenvalues(tank)
    quad   = quadrature(tank)
    τ      = pml(tank)
    dofs   = quadrature(tank).qnodes
    J_diag = Diagonal([jacobian_det(τ, dof) for dof ∈ dofs])
    S, D   = tank.S, tank.D
    # solve Ax = λBx, with A = 0.5*|J|^{-1} + D
    A = 0.5 * inv(J_diag) - D
    Γ = freesurface(tank)
    Is = dom2qtags(quad, Γ) # index of dofs on free surface
    B = zero(S)
    B[:, Is] = S[:, Is]
    @info "computing generalized eigenvalue decomposition"
    t = @elapsed begin
        F = eigen(A, B)
    end
    @info "|-- took $t seconds"
    return F
end

function solve!(tank::WaveTank, f::Function)
    quad = quadrature(tank)
    return solve!(tank, [f(dof) for dof ∈ quad.qnodes])
end

function solution(tank)
    G = tank.Sop.kernel
    dG = tank.Dop.kernel
    f = tank.f
    quad = quadrature(tank)
    dofs = quadrature(tank).qnodes
    k = impedance(tank)
    Γ = freesurface(tank)
    Is = dom2qtags(quad, Γ) # index of dofs on free surface
    τ = pml(tank)
    J_diag = Diagonal([jacobian_det(τ, dof) for dof ∈ dofs])
    f′ = J_diag * f
    σ = tank.σ
    @. f′[Is] -= k * σ[Is] # f - ω²/g |J| φ
    𝒮 = IntegralPotential(G, quad)[f′]
    𝒟 = IntegralPotential(dG, quad)[σ]
    sol = (x) -> 𝒟(x) - 𝒮(x)
    return sol
end

# exact solutions
function plane_wave(tank::WaveTank; θ = 0)
    λ = lengthscale(tank)
    d = depth(tank) / λ |> NoUnits
    β = wavenumber(tank)
    A = 1 / cosh(β * d) # normalization constant
    return _modal_solution(A, β, d)
end

function _modal_solution(A, β, d)
    ϕᵢ  = (dof) -> begin
        x = coords(dof)
        A * exp(im * β * x[1]) * cosh(β * (x[2] + d))
    end
    dϕᵢ = (dof) -> begin
        x = coords(dof)
        n = normal(dof)
        @assert norm(n) ≈ 1
        A * β * exp(im * β * x[1]) * (im * cosh(β * (x[2] + d)) * n[1] + sinh(β * (x[2] + d)) * n[2])
    end
    return ϕᵢ, dϕᵢ
end

function evanescent_wave(tank::WaveTank; n = 1, s = 1)
    @assert n ≥ 0
    d = depth(tank) / lengthscale(tank) |> NoUnits
    γ = evanescent_wavenumber(tank; n)
    A = 1 / cos(γ * d) # normalization constant
    return _modal_solution(A, s * im * γ, d)
end
