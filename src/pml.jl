abstract type AbstractPML end

"""
    OrthogonalLinearUniaxialPML

Change of variables given by `x̃ = a + c*(x-a)` for `x>a` and `x̃ = -a +
c*(x+a)` for `x<-a`. For `-a ≤ x ≤ a`, `x̃ = x`.
"""
struct OrthogonalLinearUniaxialPML <: AbstractPML
    a::Float64
    c::ComplexF64
end

function OrthogonalLinearUniaxialPML(;a,θ)
    OrthogonalLinearUniaxialPML(a,exp(im*θ))
end

function (f::OrthogonalLinearUniaxialPML)((x,y))
    a,c = f.a,f.c
    ỹ = ComplexF64(y)
    if x > a
        x̃ = a + (x-a)*c
    elseif x < -a
        x̃ = -a + (x+a)*c
    else
        x̃ = ComplexF64(x)
    end
    SVector(x̃,ỹ)
end

function jacobian_det(f::OrthogonalLinearUniaxialPML,(x,y))
    a,c = f.a,f.c
    α = abs(x) > a ? c : ComplexF64(1)
end

struct LaplacePML{N,S} <: Nystrom.AbstractPDE{N}
    complex_strecthing::S
end

LaplacePML(;dim=2,τ) = LaplacePML{dim,typeof(τ)}(τ)

getname(::LaplacePML) = "LaplacePML"

Nystrom.default_kernel_eltype(::LaplacePML)  = ComplexF64
Nystrom.default_density_eltype(::LaplacePML) = ComplexF64

complex_strecthing(op::LaplacePML) = op.complex_strecthing

function (SL::SingleLayerKernel{T,S})(target,source)::T  where {T,S<:LaplacePML}
    x = coords(target)
    y = coords(source)
    N = ambient_dimension(SL.pde)
    τ = SL.pde.complex_strecthing
    x̃ = τ(x)
    ỹ = τ(y)
    r = x̃ - ỹ
    d = sqrt(transpose(r)*r)
    d==0 && (return zero(T))
    if N==2
        return -1/(2π)*log(d)
    elseif N==3
        return 1/(4π)/d
    else
        notimplemented()
    end
end

function (DL::DoubleLayerKernel{T,S})(target,source)::T where {T,S<:LaplacePML}
    x,y,ny = coords(target), coords(source), normal(source)
    N  = ambient_dimension(DL.pde)
    τ  = DL.pde.complex_strecthing
    x̃  = τ(x)
    ỹ  = τ(y)
    r  = x̃ - ỹ
    d  = sqrt(transpose(r)*r)
    d == 0 && (return zero(T))
    if N==2
        return 1/(2π)/(d^2) * transpose(r)*ny
    elseif N==3
        return 1/(4π)/(d^3) * transpose(r)*ny
    else
        notimplemented()
    end
end


#=
    Infinite depth two-dimensional water wave Gree function
=#
struct InfiniteDepthWaterWaves{N,S} <: Nystrom.AbstractPDE{N}
    k::S # ω²/g
end

InfiniteDepthWaterWaves(;k,dim=2) = InfiniteDepthWaterWaves{dim,typeof(k)}(k)

Nystrom.default_kernel_eltype(::InfiniteDepthWaterWaves)  = ComplexF64
Nystrom.default_density_eltype(::InfiniteDepthWaterWaves) = ComplexF64

function (SL::SingleLayerKernel{T,S})(target,source)::T  where {T,S<:InfiniteDepthWaterWaves}
    N = 2
    k = SL.pde.k
    x = coords(target)
    x̄ = setindex(x,-x[N],N) # image across surface
    y = coords(source)
    r = y-x
    r̄ = y-x̄
    d = norm(r)
    d̄ = norm(r̄)
    d==0 && (return zero(T))
    if N==2
        x1,x2 = x
        y1,y2 = y
        x2,y2 = -x2,-y2
        return -1/(2π)*log(d) + 1/(2π)*log(d̄) + im*exp(-k*(y2+x2))*cos(k*(y1-x1)) - exp(-k*(y2+x2))/(2π)*(exp(im*k*(y1-x1))*SpecialFunctions.expinti(k*((y2+x2) - im*(y1-x1))) + exp(-im*k*(y1-x1))*SpecialFunctions.expinti(k*((y2+x2) + im*(y1-x1))))
    else
        notimplemented()
    end
end

function (DL::DoubleLayerKernel{T,S})(target,source)::T  where {T,S<:InfiniteDepthWaterWaves}
    N = 2
    k = DL.pde.k
    x = coords(target)
    x̄ = setindex(x,-x[N],N) # image across surface
    y = coords(source)
    r = x-y
    r̄ = x̄-y
    d = norm(r)
    d̄ = norm(r̄)
    d==0 && (return zero(T))
    ny = normal(source)
    if N==2
        x1,x2 = x
        y1,y2 = y
        x2,y2 = -x2,-y2
        v1 = y1-x1
        v2 = y2+x2
        nt = transpose(ny)
        return ((nt*r)/(2π*d^2) + (nt*r̄)/(2π*d̄^2) - (nt*SVector(sin(k*v1),-cos(k*v1)))*im*k*exp(-k*v2)
            +k/(2π)*exp(-k*v2)*(
                (nt*SVector(-im,-1))*exp(im*k*v1)*SpecialFunctions.expinti(k*(v2-im*v1)) + (nt*SVector(im,-1))*exp(-im*k*v1)*SpecialFunctions.expinti(k*(v2+im*v1))
            ))
    else
        notimplemented()
    end
end
