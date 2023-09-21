"""
    struct PhysicalParameters

Parameters with physical units. Fields include:
- `frequency`
- `gravity`
- `depth`
"""
Base.@kwdef mutable struct PhysicalParameters
    frequency::typeof(1.0u"Hz")
    gravity::typeof(1.0u"m/s^2") = 9.81 * u"m/s^2"
    depth::typeof(1.0u"m") = Inf * u"m"
end

# # mesh PhysicalParameters
# meshsize::Float64  = Inf
# qorder::Int        = 0

frequency(p::PhysicalParameters) = p.frequency
gravity(p::PhysicalParameters)   = p.gravity
depth(p::PhysicalParameters)     = p.depth

function wavelength(p::PhysicalParameters)
    ω = frequency(p)
    g = gravity(p)
    d = depth(p)
    if isinf(d)
        k₀ = ω^2 / g
    else
        k₀ = find_zero(k -> k * tanh(k * d) - ω^2 / g, ω^2 / g)
    end
    return 2π * uconvert(u"m", 1 / k₀)
end

function wavenumber(p::PhysicalParameters)
    λ = wavelength(p)
    return 2π / λ
end

function evanescent_wavelength(p::PhysicalParameters; n = 1)
    @assert n ≥ 1
    ω = frequency(p)
    g = gravity(p)
    d = depth(p)
    if isinf(d)
        @warn "evanescent wavenumber is not defined for infinite depth"
        γ = 0.0
    else
        f = γ -> γ * tan(γ * d) + ω^2 / g
        a = (1 + eps()) * (n * π - π / 2) / d
        b = (1 - eps()) * (n * π) / d
        @debug f(a), f(b)
        γ = find_zero(f, (a, b))
    end
    return 2π / γ
end

evanescent_wavenumber(p::PhysicalParameters; n = 1) = 2π / evanescent_wavelength(p; n)

# exact solutions
function plane_wave(parameters::PhysicalParameters; θ = 0)
    d = depth(parameters)
    β = wavenumber(parameters)
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

function evanescent_wave(parameters::PhysicalParameters; n = 1, s = 1)
    @assert n ≥ 0
    d = depth(parameters)
    γ = evanescent_wavenumber(parameters; n)
    A = 1 / cos(γ * d) # normalization constant
    return _modal_solution(A, s * im * γ, d)
end
