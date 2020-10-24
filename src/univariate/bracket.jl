function _bracket(f, x, step, scale, max_iter)
    a, ya = x, f(x)
    b, yb = a + step, f(a + step)
    if yb > ya
        a, b = b, a
        ya, yb = yb, ya
        step = -step
    end
    iter = 0
    while iter < max_iter
        iter += 1
        c, yc = b + step, f(b + step)
        if yc > yb
            return a < c ? (a, c) : (c, a)
        end
        a, ya, b, yb = b, yb, c, yc
        step *= scale
    end
    return nothing
end

"""
    bracket(f::Function, x::T=0.0; kwargs...) where {T<:AbstractFloat}

Create a bracket around a local minimum from the intial point `x`. Return the bracket in a 
tuple if successful, and `nothing` otherwise. 

# Keywords
- `step=0.01`: the step size to expand the bracket
- `scale=2.0`: the scale factor for step size at each iteration
- `max_iter::Integer=100`: the maximum number of iterations
"""
function bracket(f::Function, x::T=0.0; 
        step=0.01, scale=2.0, max_iter::Integer=100) where {T<:AbstractFloat}
    _bracket(f, x, T(step), T(scale), max_iter)
end