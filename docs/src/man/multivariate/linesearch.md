# Line Search Algorithms

For the mathematical background behind line search algorithms, [this Pluto notebook](https://mybinder.org/v2/gh/fonsp/pluto-on-binder/master?urlpath=pluto/open?url=https%253A%252F%252Fgist.githubusercontent.com%252Flhnguyen-vn%252F4ecbc8b9e52b82936b909771fd379e00%252Fraw%252F2cb4e7df078a016b23c0b82c005acadd0b991db5%252Flinesearch.jl) 
details a quick overview with proof of convergence. Optini's [`LineSearch`] abstraction provides
a simple interface to customize and experiment with the basic components: the search direction, 
the initial step length guess, and the method to select an appropriate step length. In particular,
a few examples from Nocedal & Wright's *Numerical Optimization* have already been implemented:
- Search direction: `Steepest`, `Newton`
- Initial guess: [`StaticInitial`](@ref), [`PreviousDecreaseInitial`](@ref), [`QuadraticInitial`](@ref)
- Step length selection method: [`StaticLineSearch`](@ref), [`ExactLineSearch`](@ref), [`BacktrackingLineSearch`](@ref), [`InterpolationLineSearch`](@ref), [`StrongWolfeLineSearch`](@ref)

## Search Direction Interface

| Function    | Input               | Output                                                                          |
|:----------- |:------------------- |:------------------------------------------------------------------------------- |
| `order`     | `YourType`          | order of information required, e.g. `Optini.SecondOrder()` for Newton direction |
| `direction` | `YourType`, `state` | the search direction                                                            |

To define a custom search direction, you must implement `order` and `direction` on your type.
The former tells Optini whether your method requires first-order and/or second-order information,
while the latter computes the direction for line search. Optini will compute each iteration's
`state` with the objective value, first-order and/or second-order information necessary for 
your algorithm. The following example, for instance, is how Optini defines the steepest 
descent direction internally:

```julia
import Optini: order, direction, FirstOrder

# Define singleton type for dispatch purposes
struct Steepest end

# Steepest descent requires first-order information
order(::Steepest) = FirstOrder()

# Since `Steepest` is now classified as `FirstOrder`, Optini will prepare `FirstOrderState`
# information at each iteration. We access the gradient in the field `∇f`, and simply define
# the search direction as the opposite of the gradient.
direction(::Steepest, state) = -state.∇f
```

## Initial Step Length Interface

| Function   | Input                         | Output                                                                           |
|:---------- |:----------------------------- |:-------------------------------------------------------------------------------- |
| `reset!`   | `YourType`                    | reset the method, especially for initial guess dependent on previous iterations  |
| `update!`  | `YourType`, `state`, `p`, `α` | update the method, especially for initial guess dependent on previous iterations |
| `YourType` | `state`, `p`                  | the initial step length guess                                                    |

A custom initial step length can be implemented by defining `reset!`, `update!`, and a functor
to compute the initial guess. For example, a fixed step length guess might be implemented as:

```julia
# For simplicity we're only defining `Float64` initial here. Internally, however, Optini 
# uses the parametric struct `StaticInitial{T}`.
struct StaticInitial
    α::Float64
end

# Static initials don't hold any information that needs to be reset. If you are subtyping 
# `Optini.AbstractInitial{T}`, this is already the fallback behavior.
function reset!(::StaticInitial) end

# Some initial methods require information about the last search direction `p` and step 
# length `α`, to choose an initial for the current iteration and therefore need to be updated. 
# `StaticInitial`, however, doesn't need to update anything here. Similar to `reset!`, this 
# is the default fallback behavior for types subtyping `Optini.AbstractInitial{T}`.
function update!(::StaticInitial, state, p, α) end

# Return the initial step length guess, given the current `state` and the direction `p`
(si::StaticInitial)(state, p) = si.α
```

## Step Length Selection Method Interface

| Function   | Input                   | Output                      |
|:---------- |:----------------------- |:--------------------------- |
| `YourType` | `f`, `state`, `p`, `α₀` | the appropriate step length |

A custom step length selection scheme receives the objective function `f`, the current `state`, 
the search direction `p`, and the initial step length `α₀` to produce the next step length. 
For instance, the following code defines a fixed line search method:

```julia
struct StaticLineSearch end

# Static line search method simply uses the initial step length
(sls::StaticLineSearch)(f, state, p, α₀) = α₀
```

## References

```@docs
LineSearch
GradientDescent
NewtonDescent
StaticInitial
PreviousDecreaseInitial
QuadraticInitial
StaticLineSearch
ExactLineSearch
BacktrackingLineSearch
InterpolationLineSearch
StrongWolfeLineSearch
```