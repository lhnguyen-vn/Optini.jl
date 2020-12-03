"""
    BacktrackingLineSearch{C, T}

`BacktrackingLineSearch` starts from an initial step length guess and gradually backtracks 
until the Armijo sufficient decrease condition is met.

# Fields
- `c::C`: the constant for the Armijo condition in the open interval (0, 1)
- `scale::T`: the scale factor to shrink the step length
"""
struct BacktrackingLineSearch{C, T}
    c::C
    scale::T

    function BacktrackingLineSearch(c, scale)
        0 < c < 1 || error("`c` must be in the open interval (0, 1)")
        0 < scale < 1 || error("Scale factor must be between 0 and 1")
        return new{typeof(c), typeof(scale)}(c, scale)
    end
end

"""
    BacktrackingLineSearch(; kwargs...)

Initiate `BacktrackingLineSearch`.

# Keywords
- `c=1e-4`: the constant for the Armijo condition in the open interval (0, 1)
- `scale=0.5`: the scale factor to shrink the step length
"""
BacktrackingLineSearch(; c=1e-4, scale=0.5) = BacktrackingLineSearch(c, scale)

function (bls::BacktrackingLineSearch)(f, state, p, α₀::T) where {T}
    x = state.x
    y = state.f
    ∇y = state.∇f
    dϕ₀ = ∇y⋅p
    c = bls.c
    scale = T(bls.scale)
    while f(x + α₀*p) > y + c*α*dϕ₀
        α₀ *= scale
    end
    return α₀
end