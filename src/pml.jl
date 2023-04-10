"""
    OrthogonalLinearPML

Change of variables given by `x̃ = a + c*(x-a)` for `x>a` and `x̃ = -a +
c*(x+a)` for `x<-a`. For `-a ≤ x ≤ a`, `x̃ = x`.
"""
struct OrthogonalLinearPML
    a::Float64
    c::ComplexF64
end

function OrthogonalLinearPML(;a,θ)
    OrthogonalLinearPML(a,exp(im*θ))
end

function (f::OrthogonalLinearPML)(dof)
    x   = WPB.coords(dof)
    N   = length(x)
    a,c = f.a,f.c
    WPB.svector(N) do d
        xd = x[d]
        if d == N || abs(xd) <= a
            Complex(xd)
        elseif xd > a
            a + (xd-a)*c
        else
            -a + (xd+a)*c
        end
    end
end

function jacobian(f::OrthogonalLinearPML,dof)
    a,c = f.a,f.c
    x = WPB.coords(dof)
    N = length(x)
    diag = svector(N) do d
        xd = x[d]
        (d == N || xd <= a) ? one(c) : c
    end
    Diagonal(diag)
end

function jacobian_det(f::OrthogonalLinearPML,dof)
    x = WPB.coords(dof)
    N = length(x)
    a,c = f.a,f.c
    prod(1:N-1) do d
        xd = x[d]
        abs(xd) > a ? c : one(c)
    end
end

struct OrthogonalQuadraticPML
    a::Float64
    c::ComplexF64
end

function OrthogonalQuadraticPML(;a,θ)
    OrthogonalQuadraticPML(a,exp(im*θ))
end

function (f::OrthogonalQuadraticPML)(dof)
    x   = WPB.coords(dof)
    N   = length(x)
    a,c = f.a,f.c
    WPB.svector(N) do d
        xd = x[d]
        if d == N || abs(xd) <= a
            Complex(xd)
        elseif xd > a
            a + (xd-a)^2*c
        else
            -a - (xd-a)^2*c
        end
    end
end

function jacobian(f::OrthogonalQuadraticPML,dof)
    a,c = f.a,f.c
    x = WPB.coords(dof)
    N = length(x)
    diag = WPB.svector(N) do d
        xd = x[d]
        (d == N || xd <= a) ? one(c) : xd > a ? 2*(xd-a)*c : -2*(xd-a)*c
    end
    Diagonal(diag)
end

function jacobian_det(f::OrthogonalQuadraticPML,dof)
    x = WPB.coords(dof)
    N = length(x)
    a,c = f.a,f.c
    prod(1:N-1) do d
        xd = x[d]
        abs(xd) <= a ? one(c) : 2*(xd-a)*c
    end
end

struct LaplacePML{N,S} <: WPB.AbstractPDE{N}
    complex_strecthing::S
end

LaplacePML(;dim=2,τ) = LaplacePML{dim,typeof(τ)}(τ)

getname(::LaplacePML) = "LaplacePML"

WPB.default_kernel_eltype(::LaplacePML)  = ComplexF64
WPB.default_density_eltype(::LaplacePML) = ComplexF64

complex_strecthing(op::LaplacePML) = op.complex_strecthing

_log(z::Complex) = isreal(z) ? log(real(z)) : log(z)

function (SL::WPB.SingleLayerKernel{T,S})(target,source)::T  where {T,S<:LaplacePML}
    x = WPB.coords(target)
    y = WPB.coords(source)
    N = WPB.ambient_dimension(SL.pde)
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
        WPB.notimplemented()
    end
end

function (DL::WPB.DoubleLayerKernel{T,S})(target,source)::T where {T,S<:LaplacePML}
    x,y,ny = WPB.coords(target), WPB.coords(source), WPB.normal(source)
    N  = WPB.ambient_dimension(DL.pde)
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
        WPB.notimplemented()
    end
end

struct LaplaceSymPML{N,S} <: WPB.AbstractPDE{N}
    complex_strecthing::S
end

LaplaceSymPML(;dim=2,τ) = LaplaceSymPML{dim,typeof(τ)}(τ)

getname(::LaplaceSymPML) = "LaplaceSymPML"

WPB.default_kernel_eltype(::LaplaceSymPML)  = ComplexF64
WPB.default_density_eltype(::LaplaceSymPML) = ComplexF64

complex_strecthing(op::LaplaceSymPML) = op.complex_strecthing

function (SL::WPB.SingleLayerKernel{T,S})(target,source)::T  where {T,S<:LaplaceSymPML}
    x = WPB.coords(target)
    y = WPB.coords(source)
    N = WPB.ambient_dimension(SL.pde)
    τ = SL.pde.complex_strecthing
    x̃ = τ(x)
    ỹ = τ(y)
    x̄ = setindex(x̃,-x̃[N],N) # image across surface
    r = x̃ - ỹ
    r̄ = x̄ - ỹ
    d = sqrt(transpose(r)*r)
    d̄ = sqrt(transpose(r̄)*r̄)
    d==0 && (return zero(T))
    if N==2
        return -1/(2π)*(log(d) - log(d̄))
    elseif N==3
        return 1/(4π)/d
    else
        WPB.notimplemented()
    end
end

function (DL::WPB.DoubleLayerKernel{T,S})(target,source)::T where {T,S<:LaplaceSymPML}
    x,y,ny = WPB>coords(target), WPB.coords(source), WPB.normal(source)
    N  = WPB.ambient_dimension(DL.pde)
    τ  = DL.pde.complex_strecthing
    x̃  = τ(x)
    ỹ  = τ(y)
    x̄ = setindex(x̃,-x̃[N],N) # image across surface
    r  = x̃ - ỹ
    d  = sqrt(transpose(r)*r)
    r̄ = x̄ - ỹ
    d̄ = sqrt(transpose(r̄)*r̄)
    d == 0 && (return zero(T))
    if N==2
        return 1/(2π)/(d^2) * transpose(r)*ny - 1/(2π)/(d̄^2) * transpose(r̄)*ny
    elseif N==3
        return 1/(4π)/(d^3) * transpose(r)*ny
    else
        WPB.notimplemented()
    end
end


#=
    Infinite depth two-dimensional water wave Gree function
=#
struct InfiniteDepthWaterWaves{N,S} <: WPB.AbstractPDE{N}
    k::S # ω²/g
end

InfiniteDepthWaterWaves(;k,dim=2) = InfiniteDepthWaterWaves{dim,typeof(k)}(k)

WPB.default_kernel_eltype(::InfiniteDepthWaterWaves)  = ComplexF64
WPB.default_density_eltype(::InfiniteDepthWaterWaves) = ComplexF64

function (SL::WPB.SingleLayerKernel{T,S})(target,source)::T  where {T,S<:InfiniteDepthWaterWaves}
    N = 2
    k = SL.pde.k
    x = WPB.coords(target)
    x̄ = setindex(x,-x[N],N) # image across surface
    y = WPB.coords(source)
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
        WPB.notimplemented()
    end
end

function (DL::WPB.DoubleLayerKernel{T,S})(target,source)::T  where {T,S<:InfiniteDepthWaterWaves}
    N = 2
    k = DL.pde.k
    x = WPB.coords(target)
    x̄ = setindex(x,-x[N],N) # image across surface
    y = WPB.coords(source)
    r = x-y
    r̄ = x̄-y
    d = norm(r)
    d̄ = norm(r̄)
    d==0 && (return zero(T))
    ny = WPB.normal(source)
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
        WPB.notimplemented()
    end
end
