function ackley(x; a=20, b=0.2, c=2π)
    d = length(x)
    return -a*exp(-b*sqrt(sum(x.^2)/d)) - 
        exp(sum(cos.(c*x))/d) + a + exp(1)
end

booth(x) = (x[1] + 2*x[2] - 7)^2 + (2*x[1] + x[2] -5)^2

function branin(x; a=1, b=5.1/(4π^2), c=5/π, r=6, s=10, t=1/(8π))
    return a*(x[2] - b*x[1]^2 + c*x[1] - r)^2 + s*(1-t)*cos(x[1]) + s
end

function flower(x; a=1, b=1, c=4)
    return a*norm(x) + b*sin(c*atan(x[2], x[1]))
end

function michalewicz(x; m=10)
    return -sum(sin(xᵢ)*sin(i*xᵢ^2/π)^(2m) for (i, xᵢ) in enumerate(x))
end

rosenbrock(x; a=1, b=100) = (a - x[1])^2 + b*(x[2] - x[1]^2)^2

wheeler(x, a=1.5) = -exp(-(x[1]*x[2] - a)^2 - (x[2] - a)^2)

function circle(x)
    θ = x[1]
    r = 0.5 + 0.5 * (2*x[2]/(1 + x[2]^2))
    return [1-r*cos(θ); 1-r*sin(θ)]
end