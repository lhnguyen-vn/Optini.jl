function optimize(f::Function, lower::T, upper::T; 
        alg::UnivariateAlgorithm=GoldenSection(), 
        rel_tol=sqrt(eps(T)), 
        abs_tol=eps(T), 
        max_iter::Integer=1_000) where {T<:AbstractFloat}
    lower < upper || error("Variable lower bound has to be smaller than upper bound.")
    _optimize(f, lower, upper, alg; rel_tol=T(rel_tol), abs_tol=T(abs_tol), max_iter)
end

"""
    optimize(f::Function, lower::Real, upper::Real; kwargs...)

Compute the minimizer of the function `f` using the bracket (`lower`, `upper`).

# Keywords
- `alg::UnivariateAlgorithm=GoldenSection()`: the algorithm choice
- `rel_tol=sqrt(eps(T))`: relative tolerance
- `abs_tol=eps(T)`: absolute tolerance
- `max_iter::Integer=1_000`: maximum numer of iterations
"""
function optimize(f::Function, lower::Real, upper::Real; kwargs...)
    T = promote_type(typeof(lower/1), typeof(upper/1))
    optimize(f, T(lower), T(upper); kwargs...)
end

