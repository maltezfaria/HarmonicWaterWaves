"""
    OrthogonalPML
"""
@kwdef struct OrthogonalPML
    pml_start::Float64
    stretch_start::Float64 = Inf
    pml_strength::Float64 = 1
    stretch_strength::Float64 = 1
end

function (f::OrthogonalPML)(dof)
    x   = coords(dof)
    N   = length(x)
    a,b,c,d = f.pml_start,f.stretch_start,f.pml_strength,f.stretch_strength
    svector(N) do dim
        xd = x[dim]
        if dim == N
            return ComplexF64(xd)
        else
            if xd > a && xd < b
                ComplexF64(xd + (xd-a)*c*im)
            elseif xd < -a && xd > -b
                ComplexF64(xd + (xd+a)*c*im)
            elseif xd >= b
                ComplexF64(d*(xd-b) + d*(xd-b)*c*im + b + (b-a)*c*im)
            elseif xd <= -b
                ComplexF64(d*(xd+b) + d*(xd+b)*c*im - b - (b-a)*c*im)
            else
                ComplexF64(xd)
            end
        end
    end
end

function jacobian_det(f::OrthogonalPML,dof)
    x = coords(dof)
    N = length(x)
    a,b,c,d = f.pml_start,f.stretch_start,f.pml_strength,f.stretch_strength
    prod(1:N-1) do dim
        xd = x[dim]
        # abs(xd) > a ? c : one(c)
        a < abs(xd) < b ? 1+im*c : abs(xd) >= b ? d + d*im*c : one(ComplexF64)
    end
end

struct LaplacePML{N,S} <: AbstractPDE{N}
    complex_strecthing::S
end

LaplacePML(;dim=2,τ) = LaplacePML{dim,typeof(τ)}(τ)

getname(::LaplacePML) = "LaplacePML"

default_kernel_eltype(::LaplacePML)  = ComplexF64
default_density_eltype(::LaplacePML) = ComplexF64

complex_strecthing(op::LaplacePML) = op.complex_strecthing

function _log(z::T)::T where {T <: Complex}
    isreal(z) ? log(real(z)) : log(z)
end

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
        return -1/(2π)*_log(d)
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
    Infinite depth two-dimensional water wave Green function
=#
struct InfiniteDepthWaterWaves{N,S} <: AbstractPDE{N}
    k::S # ω²/g
end

InfiniteDepthWaterWaves(;k,dim=2) = InfiniteDepthWaterWaves{dim,typeof(k)}(k)

default_kernel_eltype(::InfiniteDepthWaterWaves)  = ComplexF64
default_density_eltype(::InfiniteDepthWaterWaves) = ComplexF64

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
