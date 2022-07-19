"""
    mutable struct WaveTank


"""
Base.@kwdef mutable struct WaveTank
    # entities
    _freesurface::Domain    = Domain()
    pml_freesurface::Domain = Domain()
    _bottom::Domain         = Domain()
    pml_bottom::Domain      = Domain()
    obstacles::Domain       = Domain()
    # parameters
    parameters::Parameters   = Parameters()
    # discretized fields
    mesh                     = nothing
    quad                     = nothing
    Sop                      = nothing # single layer operator
    Dop                      = nothing # double layer operator
    S                        = nothing
    D                        = nothing
end


function freesurface(tank::WaveTank;include_pml=true)
    if include_pml
        return tank._freesurface ∪ tank.pml_freesurface
    else
        return tank._freesurface
    end
end

pml_freesurface(tank::WaveTank) = tank.pml_freesurface

function bottom(tank::WaveTank;include_pml=true)
    if include_pml
        return tank._bottom ∪ tank.pml_bottom
    else
        return tank._bottom
    end
end

obstacles(tank::WaveTank) = tank.obstacles

function domain(tank::WaveTank)
    freesurface(tank;include_pml=true) ∪ bottom(tank;include_pml=true) ∪ obstacles(tank)
end

parameters(t::WaveTank) = t.parameters
mesh(t::WaveTank)       = t.mesh
quadrature(t::WaveTank) = t.quad
depth(t::WaveTank)        = t |> parameters |> depth
frequency(tank::WaveTank) = tank |> parameters |> frequency
gravity(tank::WaveTank)   = tank |> parameters |> gravity
impedance(tank::WaveTank) = tank |> parameters |> impedance
wavenumber(t::WaveTank)   = t |> parameters |> wavenumber
wavelength(t::WaveTank)   = t |> parameters |> wavelength

# helper functions to generate the domain
function add_free_surface!(wavetank::WaveTank,a)
    Γf      = Domain(ParametricSurfaces.line(Point2D(a, 0), Point2D(-a, 0)))
    _add_free_surface!(wavetank,Γf,a)
    return wavetank
end
function _add_free_surface!(wavetank::WaveTank,Γf::Domain,a)
    p       = parameters(wavetank)
    p.sides = a
    isempty(wavetank._freesurface) || (@warn "replacing freesurface of wavetank")
    wavetank._freesurface = Γf
    return wavetank
end

function add_flat_bottom!(wavetank,d)
    @assert d > 0 "depth must be positive: got d=$d"
    p = parameters(wavetank)
    a = p.sides
    l = ParametricSurfaces.line(Point2D(-a, -d), Point2D(a, -d))
    Γb = Domain(l)
    add_bottom!(Γb, d, wavetank)
end

"""
    add_bottom!(wavetank,Γb,d)

Set the bottom of `wavetank` to `Γb`, which is assumed to be of constant depth
`d` as the horizontal variables `|𝐱| → ∞`.
"""
function add_bottom!(Γb::Domain, d, wavetank=WAVETANK)
    p       = parameters(wavetank)
    p.depth = d
    isempty(wavetank._bottom) || (@warn "replacing bottom of wavetank")
    wavetank._bottom = Γb
    return Γb
end

function add_obstacle!(wavetank,obs::AbstractEntity)
    push!(wavetank.obstacles, obs)
end

function add_pml_layer!(wavetank,l)
    p   = parameters(wavetank)
    a   = p.sides
    d   = p.depth
    # free surface pml
    layer_top_right = ParametricSurfaces.line(Point2D(a+l, 0), Point2D(a, 0))
    layer_top_left  = ParametricSurfaces.line(Point2D(-a, 0), Point2D(-a-l, 0))
    isempty(wavetank.pml_freesurface) || (@warn "replacing pml layer on freesurface")
    wavetank.pml_freesurface = Domain([layer_top_left,layer_top_right])
    # bottom pml (if needed)
    if d < Inf
        layer_bottom_right   = ParametricSurfaces.line(Point2D(xr+l, -d), Point2D(xr, -d))
        layer_bottom_left    = ParametricSurfaces.line(Point2D(xl, -d), Point2D(xl-l, -d))
        isempty(wavetank.pml_bottom) || (@warn "replacing pml layer on bottom")
        wavetank.pml_bottom  = Domain([layer_bottom_left,layer_bottom_right])
    end
    return wavetank
end

function discretize!(tank::WaveTank;meshsize,order)
    Γ    = domain(tank)
    @info "meshing wavetank"
    mesh = ParametricSurfaces.meshgen(Γ;meshsize)
    tank.mesh = mesh
    @info "creating a quadrature for the mesh"
    quad   = NystromMesh(mesh;order)
    tank.quad = quad
    # assemble lazy integral operators
    @info "creating lazy integral operators"
    p       = parameters(tank)
    a       = p.sides
    θ       = p.θ
    τ       = OrthogonalLinearUniaxialPML(;a,θ)
    op      = LaplacePML(;dim=2,τ)
    Sop     = SingleLayerOperator(op,quad)
    Dop     = DoubleLayerOperator(op,quad)
    tank.Sop = Sop
    tank.Dop = Dop
    return tank
end

function assemble_operators!(tank::WaveTank;correction=:quadgk)
    m,n = size(tank.Sop)
    @info "assembling $m × $n single- and double-layer operators"
    if correction == :quadgk
        tank.S       = tank.Sop |> Nystrom.assemble_gk
        tank.D       = tank.Dop |> Nystrom.assemble_gk
    elseif correction == :none
        tank.S       = tank.Sop |> Matrix
        tank.D       = tank.Dop |> Matrix
    else
        error("unrecognized `correction` option")
    end
    return tank
end

function solve(tank::WaveTank,f::AbstractVector)
    p       = parameters(tank)
    quad    = quadrature(tank)
    a       = p.sides
    θ       = p.θ
    dofs    = quadrature(tank).dofs
    k       = impedance(tank)
    J_diag  = Diagonal([abs(dof.coords[1]) < a ? one(ComplexF64) : exp(im*θ) for dof in dofs])
    @assert J_diag*f ≈ f "source term must vanish inside the PML"
    # compute L = 0.5*|J| + D - ω^2/g * Sfs
    S,D = tank.S, tank.D
    L   = 0.5*inv(J_diag) + D
    Γ   = freesurface(tank;include_pml=true)
    Is  = Nystrom.dom2dof(quad,Γ) # index of dofs on free surface
    @. L[:,Is] = L[:,Is] - k*S[:,Is]
    rhs = S*f
    σ   = L\rhs
    ϕ   = inv(J_diag)*σ
    return Density(ϕ,quad)
end

function solve(tank::WaveTank,f::Function)
    quad = quadrature(tank)
    solve(tank,[f(dof) for dof in quad.dofs])
end

# plotting

@recipe function f(tank::WaveTank; meshsize=0.1)
    @series begin
        lc := :blue
        label := ""
        meshsize := meshsize
        tank._freesurface
    end
    @series begin
        lc := :lightblue
        label := ""
        meshsize := meshsize
        tank.pml_freesurface
    end
    @series begin
        lc := :black
        label := "bottom"
        meshsize := meshsize
        tank._bottom
    end
    @series begin
        label := "obstacles"
        lc := :red
        meshsize := meshsize
        tank.obstacles
    end
end
