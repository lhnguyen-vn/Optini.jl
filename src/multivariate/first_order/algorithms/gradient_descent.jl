"""
    GradientDescent{T<:AbstractLineSearch} <: FirstOrderAlgorithm

`GradientDescent` uses steepest descent and line search to minimize the objective function.

# Fields
- `linesearch::T`: line search method
"""
struct GradientDescent{T<:AbstractLineSearch} <: FirstOrderAlgorithm
    linesearch::T
end

"""
    GradientDescent(; linesearch=StaticLineSearch())

Initiate `GradientDescent` algorithm.
"""
GradientDescent(; linesearch=StaticLineSearch()) = GradientDescent(linesearch)

(gd::GradientDescent)(state) = -state.∇f