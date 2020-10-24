"""
    GradientDescent{T<:AbstractLineSearch} <: FirstOrderAlgorithm

`GradientDescent` uses steepest descent and line search to minimize the objective function.

# Fields
-`linesearch::T`: line search method
"""
struct GradientDescent{T<:AbstractLineSearch} <: FirstOrderAlgorithm
    linesearch::T
    
    @doc """
        GradientDescent([line_search=StaticLineSearch()])
    
    Initiate `GradientDescent` algorithm.
    """
    GradientDescent(linesearch=StaticLineSearch()) = new{typeof(linesearch)}(linesearch)
end

(gd::GradientDescent)(state) = -state.∇f