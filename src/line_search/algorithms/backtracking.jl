"""
    struct BacktrackingLineSearch{I<:AbstractInitial, C, T} <: AbstractLineSearch

`BacktrackingLineSearch` starts from an initial step length guess and gradually backtracks 
until the Armijo sufficient decrease condition is met.

# Fields
- `init::I`: the initial step length method
- `c::C`: the constant for the Armijo condition in the open interval (0, 1)
- `scale::T`: the scale factor to shrink the step length
"""
struct BacktrackingLineSearch{I<:AbstractInitial, C, T} <: AbstractLineSearch
    init::I
    c::C
    scale::T

    function BacktrackingLineSearch(init, c, scale)
        0 < c < 1 || error("`c` must be in the open interval (0, 1)")
        return new{typeof(init), typeof(c), typeof(scale)}(init, c, scale)
    end
end

"""
    BacktrackingLineSearch(; kwargs...)

Initiate `BacktrackingLineSearch`.

# Keywords
- `init=StaticInitial()`: the initial step length method
- `c=1e-4`: the constant for the Armijo condition in the open interval (0, 1)
- `scale=0.5`: the scale factor to shrink the step length
"""
function BacktrackingLineSearch(;
        init=StaticInitial(), 
        c=1e-4,
        scale=0.5)
    BacktrackingLineSearch(init, c, scale)
end

function (bls::BacktrackingLineSearch{<:AbstractInitial{T}})(f, state, p) where {T}
    x = state.x
    y = state.f
    ∇y = state.∇f
    dϕ₀ = ∇y⋅p
    α = bls.init(state, p)
    c = bls.c
    scale = T(bls.scale)
    while f(x + α*p) > y + c*α*dϕ₀
        α *= scale
    end
    return α
end