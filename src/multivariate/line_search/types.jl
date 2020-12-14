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

"""
    LineSearch(; kwargs...)

Initiate `LineSearch` algorithm.

# Keywords
- `direction=Steepest()`: determines the search direction
- `initial=StaticInitial(1.0)`: the method to guess an initial step length
- `method=BacktrackingLineSearch()`: the method to compute the appropriate next step
"""
function LineSearch(; 
        direction=Steepest(),
        initial=StaticInitial(),
        method=BacktrackingLineSearch())
    LineSearch(direction, initial, method)
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
- `initial=StaticInitial(0.001)`: starts line search with an initial step length guess
- `method=BacktrackingLineSearch()`: computes the appropriate next step
"""
function GradientDescent(;
        initial=StaticInitial(0.001), 
        method=BacktrackingLineSearch())
    GradientDescent(initial, method)
end

"""
    GradientDescent(initial, method)

Initiate `GradientDescent`.

# Arguments
- `initial`: starts line search with an initial step length guess
- `method`: computes the appropriate next step
"""
GradientDescent(initial, method) = LineSearch(Steepest(), initial, method)

"""
    NewtonDescent{I, M}

Aliasing `LineSearch` with Newton method.
"""
const NewtonDescent{I, M} = LineSearch{Newton, I, M} where {I, M}

"""
    NewtonDescent(; kwargs...)

Initiate `NewtonDescent`.

# Keywords
- `initial=QuadraticInitial(0.001)`: starts line search with an initial step length guess
- `method=BacktrackingLineSearch()`: computes the appropriate next step
"""
function NewtonDescent(;
        initial=QuadraticInitial(0.001), 
        method=BacktrackingLineSearch())
    NewtonDescent(initial, method)
end

"""
    NewtonDescent(initial, method)

Initiate `NewtonDescent`.

# Arguments
- `initial`: starts line search with an initial step length guess
- `method`: computes the appropriate next step
"""
NewtonDescent(initial, method) = LineSearch(Newton(), initial, method)