struct QuadraticFit <: UnivariateAlgorithm end

function _optimize(f, lower, mid, upper, algorithm::QuadraticFit; rel_tol, abs_tol, max_iter, max_eval)
    y_lower, y_upper, y_mid = f(lower), f(upper), f(mid)
    iter = 0
    eval = 3
    while (iter < max_iter && eval < max_eval)
        iter += 1
        eval += 1
        x = 0.5 * (y_lower * (mid^2 - upper^2) +
                y_mid * (upper^2 - lower^2) + 
                y_upper * (lower^2 - mid^2)) /
            (y_lower * (mid - upper) + y_mid * (upper - lower) + y_upper * (lower - mid))
        yx = f(x)
        x_tol = rel_tol * abs(x) + abs_tol
        mid_point = (lower + upper) / 2
        if abs(x - mid_point) <= 2x_tol - abs(upper - lower) / 2
            return UnivariateSolution(true, iter, x, yx)
        end
        iter == max_iter && return UnivariateSolution(false, iter, x, yx)
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
end

function _optimize(f, lower::T, upper::T, algorithm::QuadraticFit; 
        mid=(lower + upper) / 2, 
        rel_tol, 
        abs_tol, 
        max_iter, 
        max_eval::Integer=max_iter+3, 
        kwargs...) where {T<:AbstractFloat}
    _optimize(f, lower, T(mid), upper, algorithm; rel_tol, abs_tol, max_iter, max_eval, kwargs...)
end