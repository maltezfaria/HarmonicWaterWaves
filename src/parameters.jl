Base.@kwdef mutable struct Parameters
    frequency::Float64 = sqrt(π) # frequency
    gravity::Float64   = 1  #
    θ::Float64         = π / 4  # pml parameter
    depth::Float64     = Inf
    sides::Float64     = 0
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

wavelength(p::Parameters) = 2π/wavenumber(p)
