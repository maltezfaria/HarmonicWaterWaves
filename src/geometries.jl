#=
    Some simple geometries for testing purposes. Mostly to centralize and
    simplify the scripts in scripts_final, but may be useful as a starting point
    to build more complex shapes.
=#

"""
    wavemaker_geometry!(tank::WaveTank,l, a=0)

Create a wavemaker geometry given by `[0,l] × [-d,0]` and append it to `tank`.
The paramter `0 ≤ a < 1` controls the amplitude, relative to `d`, of a
sinusoidal *bump* on the bottom.

Special care is taken to ensure that the geometrical entities respect the
various discontinuities that may appear in the `pml_func` of `tank` (which must
be initialized beforehand).
"""
function add_wavemaker!(tank::WaveTank, l, depth_func)
    d = depth(tank) / lengthscale(tank) |> NoUnits # dimensionless depth
    τ = pml(tank)
    isnothing(τ) && error("`pml_func` must be initialized before calling this function")
    topography = WPB.ParametricEntity(0, τ.pml_start) do (s,)
        t = τ.pml_start - s
        return SVector(t, depth_func(t))
    end
    skip = isapprox(τ.pml_start, 0) # do you need a region without PML?
    if l <= τ.pml_start
        @error "wavemaker geometry truncated before the start of PML layer"
    elseif l <= τ.stretch_start
        skip || add_freesurface!(tank, 0, τ.pml_start)
        add_freesurface!(tank, τ.pml_start, l)
        skip || add_bottom!(tank, Domain(topography))
        # add_bottom!(tank,0,τ.pml_start)
        add_bottom!(tank, τ.pml_start, l)
    else
        skip || add_freesurface!(tank, 0, τ.pml_start)
        add_freesurface!(tank, τ.pml_start, τ.stretch_start)
        add_freesurface!(tank, τ.stretch_start, l)
        skip || add_bottom!(tank, Domain(topography))
        # add_bottom!(tank,0,τ.pml_start)
        add_bottom!(tank, τ.pml_start, τ.stretch_start)
        add_bottom!(tank, τ.stretch_start, l)
    end
    side = WPB.line((0, -d), (0, 0))
    add_obstacles!(tank, WPB.Domain(side))
    return tank
end
function add_wavemaker!(tank::WaveTank, l)
    d = depth(tank) / lengthscale(tank) |> NoUnits # dimensionless depth
    depth_func = (t) -> -d
    return add_wavemaker!(tank, l, depth_func)
end

function jellyfish(r₀, x₀, θ)
    rot = @SMatrix [cos(θ) -sin(θ); sin(θ) cos(θ)]
    WPB.ParametricEntity(0, 2π) do (s,)
        r = r₀ * (1 + 0.3 * cos(4 * (s + 0.5 * sin(s))))
        pt = rot * (SVector(cos(s) * r, sin(s) * r)) .+ x₀
        return pt
    end
end

function add_waveguide!(tank::WaveTank, l, depth_func = (s) -> -depth(tank))
    τ = pml(tank)
    isnothing(τ) && error("`pml_func` must be initialized before calling this function")
    a = τ.pml_start
    add_freesurface!(tank, -a, a)
    add_freesurface!(tank, a, l)
    add_freesurface!(tank, -l, -a)

    # depth between [-a,a]
    topo = WPB.ParametricEntity(-a, a) do (s,)
        # flip orientation so that normal points up, so s → -s
        return SVector(-s, depth_func(-s))
    end |> WPB.Domain
    add_bottom!(tank, topo)
    add_bottom!(tank, a, l)
    add_bottom!(tank, -l, -a)
    d = depth(tank)
    return tank
end
