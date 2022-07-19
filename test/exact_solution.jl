####Analytical solution####
##INFINITE DEPTH#####

using QuadGK

function solution_sup1(x::Real,y::Real,k::Real)
    residu(x::Real,y::Real,k::Real) = ffourier(k) * exp(k * y) * exp(im * x * k)
    g(x::Real,y::Real,t::Real,k::Real) = -im / (sqrt(2 * pi) * t^2) * (2 * exp(-x * t) - exp(-t * (x + 1)) - exp(-t * (x - 1)) ) * ( exp(im * t * y) / (im * t - k) + exp(-im * t * y) / (im * t + k) )
    ffourier(x::Real) = ( 2 - (exp(im * x) + exp(-im * x)) ) / (sqrt(2 * pi) * x^2)

    (2 * im * pi * residu(x, y, k) + quadgk(t -> g(x, y, t, k), 0, Inf)[1]) / sqrt(2 * pi)
end

function solution_inf1(x::Real, y::Real, k::Real)
    ############X<1###############
    sdpi=1/sqrt(2*pi)

    resg1=(x,y)-> sdpi*(2*exp(im*k*x)-exp(im*k*(x+1)))*exp(k*y)/(k^2)
    resg2=(x,y)-> sdpi*exp(-im*k*(x-1))*exp(k*y)/(k^2)

    f1=(t,x,y)-> im*sdpi*( 2*exp(im*x*(im*t+k/2)) - exp(im*(x+1)*(im*t+k/2)) )*exp((im*t+k/2)*y)/( (im*t-k/2)*(im*t+k/2)^2 )
    f2=(t,x,y)-> im*sdpi*( 2*exp(im*x*(im*t-k/2)) - exp(im*(x+1)*(im*t-k/2)) )*exp((-im*t+k/2)*y)/( (im*t+k/2)*(im*t-k/2)^2 )
    f3=(t,x,y)-> -im*sdpi*exp(im*(x-1)*(im*t+k/2))*exp((im*t+k/2)*y)/( (im*t-k/2)*(im*t+k/2)^2 )
    f4=(t,x,y)-> -im*sdpi*exp(im*(x-1)*(im*t-k/2))*exp((-im*t+k/2)*y)/( (im*t+k/2)*(im*t-k/2)^2 )
    fm=(t,x,y)-> sdpi*(2*exp(im*t*x)-exp(im*t*(x+1))-exp(im*t*(x-1)))*exp(sqrt(t^2)*y)/((sqrt(t^2)-k)*t^2)

    # inte1= (x,y)-> integral[@(t) f1(t,x,y)+f2(t,x,y),0,Inf,"ArrayValued", true];
    inte1= (x,y)-> quadgk(t->f1(t,x,y)+f2(t,x,y),0,Inf)[1]
    inte2= (x,y)-> quadgk(t->f3(t,x,y)+f4(t,x,y),0,-Inf)[1]
    intem= (x,y)-> quadgk(t->fm(t,x,y),-k/2,0,k/2)[1]

    ss=sdpi*real(2*im*pi*(resg1(x,y)-resg2(x,y))+inte1(x,y)+inte2(x,y)+intem(x,y))
    return ss

end

function solution_totale(x::Real, y::Real, k::Real)
    if 0 ≤ x ≤ 1
        return real(solution_inf1(x, y, k))
    elseif x > 1
        return real(solution_sup1(x, y, k))
    end
    return solution_totale(-x, y, k)
end
