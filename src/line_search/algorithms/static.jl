"""
    StaticLineSearch{T<:AbstractInitial}

`StaticLineSearch` directly uses the initial guess as the new step length.

# Fields
- `init::T`: the initial step length method
"""
struct StaticLineSearch{T<:AbstractInitial} <: AbstractLineSearch
    init::T
end

"""
    StaticLineSearch(; init=StaticInitial())

Initiates `StaticLineSearch`, initial guess defaults to `StaticInitial`.
"""
StaticLineSearch(; init=StaticInitial()) = StaticLineSearch{typeof(init)}(init)

(sls::StaticLineSearch)(f, state, p) = sls.init(state, p)