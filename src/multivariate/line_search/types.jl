abstract type AbstractInitial{T} end

"""
    LineSearch{D, I, M} <: MultivariateAlgorithm

`LineSearch` optimizes the objective function by first choosing a direction, then search 
along this direction for an appropriate step. 

# Fields
- `direction::D`: determines the search direction
- `initial::I`: the method to guess an initial step length
- `method::`: the method to compute the appropriate next step
"""
struct LineSearch{D, I, M} <: MultivariateAlgorithm
    direction::D
    initial::I
    method::M
end

order(ls::LineSearch) = order(ls.direction)

"""
    GradientDescent{I, M}

Aliasing `LineSearch` with steepest descent.
"""
const GradientDescent{I, M} = LineSearch{Steepest, I, M} where {I, M}

"""
    GradientDescent(; kwargs...)

Initiate `GradientDescent`.

# Keywords
- `initial`: starts line search with an initial step length guess
- `method`: computes the appropriate next step
"""
function GradientDescent(;
        initial=StaticInitial(0.001), 
        method=StaticLineSearch())
    LineSearch(Steepest(), initial, method)
end