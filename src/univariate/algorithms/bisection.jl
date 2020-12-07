"""
    Bisection <: UnivariateAlgorithm

Bisection search is a root-finding algorithm by maintainining a bracket of opposing 
derivative signs. At each iteration the algorithm shrinks the bracket by half until a
critical point is located or the bracket is sufficiently small. 
"""
struct Bisection <: UnivariateAlgorithm end

function _optimize(f, g, lower::T, upper::T, alg::Bisection; 
        reltol, abstol, maxiter) where {T}
    y_lower, y_upper = g(lower), g(upper)
    sign(y_lower) == sign(y_upper) && 
        error("Bisection requires initial bracket with different derivative signs.")
    y_lower == 0 && (upper = lower)
    y_upper == 0 && (lower = upper)
    x = T(NaN)
    converged = false
    iter = 0
    while iter < maxiter
        iter += 1
        x = (lower + upper) / 2
        x_tol = reltol * abs(x) + abstol
        y = g(x)
        if upper - lower < 2x_tol || y == 0
            converged = true
            break
        end
        if sign(y) == sign(y_lower)
            lower = x
        else
            upper = x
        end
    end
    return Solution(converged, iter, x, f(x))
end

function _optimize(f, lower, upper, alg::Bisection; kwargs...)
    g(x) = Zygote.gradient(f, x)[1]
    _optimize(f, g, lower, upper, alg; kwargs...)
end