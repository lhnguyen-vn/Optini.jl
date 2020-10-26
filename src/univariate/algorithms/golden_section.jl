"""
    GoldenSection <: UnivariateAlgorithm

The golden section search algorithm optimize a univariate function by iteratively shrinking 
the initial bracket around the minimum with the golden ratio as an approximation to the 
Fibonacci sequence.
"""
struct GoldenSection <: UnivariateAlgorithm end

function _optimize(f, lower::T, upper::T, alg::GoldenSection; 
        rel_tol, abs_tol, max_iter) where {T}
    p = φ - 1
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
        new_x = T(p * upper + (1 - p) * lower)
        new_yx = f(new_x)
        if yx > new_yx
            lower, x, yx = x, new_x, new_yx
        else
            lower, upper = new_x, lower
        end
    end
    return Solution(converged, iter, x, yx)
end