"""
    optimize(f::Function, lower::T, upper::T; kwargs...) where {T<:AbstractFloat}

Compute the minimizer of the function `f` using the bracket (`lower`, `upper`) and return a
`Solution` struct.

# Keyword Arguments
- `algorithm::BracketingAlgorithm=GoldenSection()`: the bracketing algorithm choice
- `rel_tol`: relative tolerance
- `abs_tol`: absolute tolerance
- `max_ter`: maximum numer of iterations
"""
function optimize(f::Function, lower::T, upper::T; 
        algorithm::BracketingAlgorithm=GoldenSection(), 
        rel_tol=sqrt(eps(T)), 
        abs_tol=eps(T), 
        max_iter::Integer=100, 
        kwargs...) where {T<:AbstractFloat}
    lower < upper || 
        throw(ErrorException("Variable lower bound has to be lower than upper bound."))
    _optimize(f, lower, upper, algorithm; 
        rel_tol=T(rel_tol), abs_tol=T(abs_tol), max_iter, kwargs...)
end

"""
    optimize(f::Function, initial::T=0.0; kwargs...) where {T<:AbstractFloat}

Compute the minimizer of the function `f` given an initial value and return a `Solution`
struct. 

# Keyword Arguments
- `algorithm::BracketingAlgorithm=GoldenSection()`: the algorithm choice
"""
function optimize(f::Function, initial::T=0.0; 
        algorithm::BracketingAlgorithm=GoldenSection(), 
        kwargs...) where {T<:AbstractFloat}
    initial_bracket = bracket(f, initial; kwargs...)
    isnothing(initial_bracket) && return BracketingSolution(false, 0, nothing, nothing)
    lower, upper = initial_bracket
    optimize(f, lower, upper; algorithm, kwargs...)
end
