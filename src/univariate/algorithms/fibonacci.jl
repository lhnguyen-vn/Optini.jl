struct Fibonacci <: UnivariateAlgorithm end

function _optimize(f, lower, upper, algorithm::Fibonacci, ϵ; rel_tol, abs_tol, max_iter)
    s = (1 - √5) / (1 + √5)
    p = 1 / (φ * (1 - s^(max_iter + 1)) / (1 - s^max_iter))
    x = p * lower + (1 - p) * upper
    yx = f(x)
    converged = false
    iter = 1
    while iter < max_iter
        x_tol = rel_tol * abs(x) + abs_tol
        mid_point = (lower + upper) / 2
        if abs(x - mid_point) <= 2x_tol - abs(upper - lower) / 2
            converged = true
            break
        end
        iter += 1
        if iter == max_iter
            new_x = ϵ * upper + (1 - ϵ) * x
        else
            new_x = p * upper + (1 - p) * lower
        end
        new_yx = f(new_x)
        if yx > new_yx
            lower, x, yx = x, new_x, new_yx
        else
            lower, upper = new_x, lower
        end
        p = 1 / (φ * (1 - s^(max_iter - iter + 2)) / (1 - s^(max_iter - iter + 1)))
    end
    return UnivariateSolution(converged, iter, x, yx)
end

function _optimize(f, lower::T, upper::T, algorithm::Fibonacci; 
        ϵ=0.01, rel_tol, abs_tol, max_iter, kwargs...) where {T<:AbstractFloat}
    _optimize(f, lower, upper, algorithm, T(ϵ); rel_tol, abs_tol, max_iter)
end