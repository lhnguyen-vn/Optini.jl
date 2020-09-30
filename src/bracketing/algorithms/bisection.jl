struct Bisection <: BracketingAlgorithm end

function _optimize(f, lower::T, upper::T, algorithm::Bisection; 
        rel_tol, abs_tol, max_iter, kwargs...) where {T<:AbstractFloat}
    f′(x) = Zygote.gradient(f, x)[1]
    y_lower, y_upper = f′(lower), f′(upper)
    sign(y_lower) == sign(y_upper) && 
        throw(ErrorException("Bisection requires initial bracket with different derivative signs."))
    y_lower == 0 && (upper = lower)
    y_upper == 0 && (lower = upper)
    iter = 1
    while iter < max_iter
        x = (lower + upper) / 2
        x_tol = rel_tol * abs(x) + abs_tol
        y = f′(x)
        if upper - lower <= x_tol || y == 0
            return BracketingSolution(true, iter, x, f(x))
        end
        iter == max_iter && return BracketingSolution(false, iter, x, f(x))
        iter += 1
        if sign(y) == sign(y_lower)
            lower = x
        else
            upper = x
        end
    end
end