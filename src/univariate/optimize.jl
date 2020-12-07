function optimize(f::Function, lower::T, upper::T; 
        alg::UnivariateAlgorithm=GoldenSection(), 
        reltol=sqrt(eps(T)), 
        abstol=eps(T), 
        maxiter::Integer=1_000) where {T<:AbstractFloat}
    lower < upper || error("Variable lower bound has to be smaller than upper bound.")
    _optimize(f, lower, upper, alg; reltol=T(reltol), abstol=T(abstol), maxiter)
end

"""
    optimize(f::Function, lower::Real, upper::Real; kwargs...)

Compute the minimizer of the function `f` using the bracket (`lower`, `upper`).

# Keywords
- `alg::UnivariateAlgorithm=GoldenSection()`: the algorithm choice
- `reltol=sqrt(eps(T))`: relative tolerance
- `abstol=eps(T)`: absolute tolerance
- `maxiter::Integer=1_000`: maximum numer of iterations
"""
function optimize(f::Function, lower::Real, upper::Real; kwargs...)
    T = promote_type(typeof(lower/1), typeof(upper/1))
    optimize(f, T(lower), T(upper); kwargs...)
end