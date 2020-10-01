"""
    GoldenSection <: BracketingAlgorithm

The golden section search algorithm optimize a univariate function by iteratively shrinking 
the initial bracket around the minimum with the golden ratio as an approximation to the 
Fibonacci sequence.
"""
struct GoldenSection <: BracketingAlgorithm end

function _optimize(f, lower, upper, algorithm::GoldenSection; rel_tol, abs_tol, max_iter, kwargs...)
    p = φ - 1
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
        new_x = p * upper + (1 - p) * lower
        new_yx = f(new_x)
        if yx > new_yx
            lower, x, yx = x, new_x, new_yx
        else
            lower, upper = new_x, lower
        end
    end
    return Solution(converged, iter, x, yx)
end