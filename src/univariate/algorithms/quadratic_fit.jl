"""
    QuadraticFit{T} <: UnivariateAlgorithm

The quadratic fit algorithm constructs a quadratic function from three points at each
iteration as an approximation to the objective. The analytical solution to the quadratic 
model is examined as the new minimizer candidate to update the three points. 
"""
struct QuadraticFit{T} <: UnivariateAlgorithm 
    ϵ::T

    function QuadraticFit{T}(ϵ) where {T}
        0 < ϵ < 1 ? new{T}(ϵ) : error("`ϵ` must be in the interval (0, 1)")
    end
end

"""
    QuadraticFit([ϵ=0.01])

Initialize quadratic fit search algorithm with optional argument `ϵ` in the interval 
(0, 1) to choose a midpoint as `mid = ϵ * lower + (1 - ϵ) * upper`.
"""
QuadraticFit(ϵ=0.5) = QuadraticFit{typeof(ϵ)}(ϵ)

function _optimize(f, lower::T, upper::T, alg::QuadraticFit; 
        reltol, abstol, maxiter) where {T}
    ϵ = alg.ϵ
    mid = T(ϵ * lower + (1 - ϵ) * upper)
    y_lower, y_upper, y_mid = f(lower), f(upper), f(mid)
    converged = false
    iter = 0
    x = T(NaN)
    yx = (typeof(y_lower))(NaN)
    while iter < maxiter
        iter += 1
        x = T(0.5 * (y_lower * (mid^2 - upper^2) + y_mid * (upper^2 - lower^2) + 
            y_upper * (lower^2 - mid^2)) /
            (y_lower * (mid - upper) + y_mid * (upper - lower) + y_upper * (lower - mid)))
        yx = f(x)
        x_tol = reltol * abs(x) + abstol
        if max(upper - x, x - lower) < 2x_tol
            converged = true
            break
        end
        if x > mid
            if yx > y_mid
                upper, y_upper = x, yx
            else
                lower, y_lower, mid, y_mid = mid, y_mid, x, yx
            end
        elseif x < mid
            if yx > y_mid
                lower, y_lower = x, yx
            else
                upper, y_upper, mid, y_mid = mid, y_mid, x, yx
            end
        end
    end
    return Solution(converged, iter, x, yx)
end