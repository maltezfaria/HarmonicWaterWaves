Base.@kwdef mutable struct Parameters
    frequency::Float64 = sqrt(π) # frequency
    gravity::Float64   = 1  #
    depth::Float64     = Inf
    # mesh parameters
    meshsize::Float64  = Inf
    qorder::Int        = 0
end

frequency(p::Parameters) = p.frequency
gravity(p::Parameters) = p.gravity
impedance(p::Parameters) = frequency(p)^2 / gravity(p)
depth(p::Parameters) = p.depth

function wavenumber(p::Parameters)
    ω = frequency(p)
    g = gravity(p)
    d = depth(p)
    k₀ = ω^2/g
    if d == Inf
        return k₀
    else
        return find_zero(k -> k * tanh(k * d) - ω^2 / g, ω^2 / g)
    end
end

function evanescent_wavenumber(p::Parameters;n=1)
    @assert n ≥ 1
    ω = frequency(p)
    g = gravity(p)
    d = depth(p)
    if d == Inf
        @warn "evanescent wavenumber is not defined for infinite depth"
        γ = 0.0
    else
        f = γ -> γ * tan(γ*d) + ω^2 / g
        a = (n*π - π/2)/d + 1e-10
        b = (n*π + π/2)/d - 1e-10
        @debug f(a),f(b)
        γ = find_zero(f,(a,b))
    end
    return γ
end

wavelength(p::Parameters) = 2π/wavenumber(p)
