function optimize(f::Function, lower::T, upper::T; 
        algorithm::BracketingAlgorithm, 
        rel_tol=sqrt(eps(T)), 
        abs_tol=eps(T), 
        max_iter::Integer=100, 
        kwargs...) where {T<:AbstractFloat}
    lower < upper || 
        throw(ErrorException("Variable lower bound has to be lower than upper bound."))
    _optimize(f, lower, upper, algorithm; 
        rel_tol=T(rel_tol), abs_tol=T(abs_tol), max_iter, kwargs...)
end

function optimize(f::Function, initial; algorithm::BracketingAlgorithm, kwargs...)
    initial_bracket = bracket(f, initial; kwargs...)
    isnothing(initial_bracket) && return BracketingSolution(false, 0, nothing, nothing)
    lower, upper = initial_bracket
    optimize(f, lower, upper; algorithm, kwargs...)
end
