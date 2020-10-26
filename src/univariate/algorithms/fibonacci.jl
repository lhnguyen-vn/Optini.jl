"""
    Fibonacci{T} <: UnivariateAlgorithm

The Fibonacci search algorithm optimize a univariate function by iteratively shrinking the 
initial bracket around the minimum. The ratio is computed from the Fibonacci sequence.
"""
struct Fibonacci{T} <: UnivariateAlgorithm 
    ϵ::T

    function Fibonacci{T}(ϵ) where {T}
        0 < ϵ < 1 ? new{T}(ϵ) : error("`ϵ` must be in the interval (0, 1)")
    end
end

"""
    Fibonacci([ϵ=0.01])

Initialize Fibonacci search algorithm with optional argument `ϵ` in the interval (0, 1)
to control how close the last iterate is to the previous one. Refer to Kochenderfer and 
Wheeler's "Algorithms for Optimization" for more information.
"""
Fibonacci(ϵ=0.01) = Fibonacci{typeof(ϵ)}(ϵ)

function _optimize(f, lower::T, upper::T, alg::Fibonacci, ϵ; 
        rel_tol, abs_tol, max_iter) where {T}
    s = (1 - √5) / (1 + √5)
    p = 1 / (φ * (1 - s^(max_iter + 1)) / (1 - s^max_iter))
    x = T(p * lower + (1 - p) * upper)
    yx = f(x)
    converged = false
    iter = 1
    while true
        x_tol = rel_tol * abs(x) + abs_tol
        if abs(upper - lower) < 2x_tol
            converged = true
            break
        end
        iter == max_iter && break
        iter += 1
        if iter == max_iter
            new_x = T(ϵ * upper + (1 - ϵ) * x)
        else
            new_x = T(p * upper + (1 - p) * lower)
        end
        new_yx = f(new_x)
        if yx > new_yx
            lower, x, yx = x, new_x, new_yx
        else
            lower, upper = new_x, lower
        end
        p = 1 / (φ * (1 - s^(max_iter - iter + 2)) / (1 - s^(max_iter - iter + 1)))
    end
    return Solution(converged, iter, x, yx)
end

function _optimize(f, lower::T, upper::T, alg::Fibonacci; 
        rel_tol, abs_tol, max_iter) where {T}
    _optimize(f, lower, upper, alg, T(alg.ϵ); rel_tol, abs_tol, max_iter)
end