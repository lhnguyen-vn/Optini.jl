"""
    InterpolationLineSearch{C, T}

`InterpolationLineSearch` uses quadratic and cubic interpolations to compute the step length
satisfying the Armijo sufficient decrease condition.

# Fields
- `c::C`: the constant for the Armijo condition in the open interval (0, 1)
- `ϵ::T`: the minimum step length decrease from the initial step length
"""
struct InterpolationLineSearch{C, T}
    c::C
    ϵ::T

    function InterpolationLineSearch(c, ϵ)
        0 < c < 1 || error("`c` must be in the open interval (0, 1)")
        return new{typeof(c), typeof(ϵ)}(c, ϵ)
    end
end

"""
    InterpolationLineSearch(; kwargs...)

Initiate `InterpolationLineSearch`.

# Keywords
- `c=1e-4`: the constant for the Armijo condition in the open interval (0, 1)
- `ϵ=1e-6`: the minimum step length decrease from the initial step length
"""
InterpolationLineSearch(; c=1e-4, epsilon=1e-6) = InterpolationLineSearch(c, epsilon)

function (ils::InterpolationLineSearch)(f, state, p, α₀::T) where {T}
    x = state.x
    ϕ₀ = state.f
    ∇y = state.∇f
    dϕ₀ = ∇y⋅p
    αₖ = α₀
    c = ils.c
    ϵ = ils.ϵ
    ϕ = α -> f(x + α * p)
    ϕαₖ = ϕ(αₖ)
    αₖ₊₁ = T(-(dϕ₀ * αₖ^2) / (2 * (ϕαₖ - ϕ₀ - dϕ₀ * αₖ)))
    while f(x + αₖ₊₁*p) > ϕ₀ + c*αₖ₊₁*dϕ₀
        ϕαₖ₊₁ = ϕ(αₖ₊₁)
        a, b = [αₖ^2 -αₖ₊₁^2; -αₖ^3 αₖ₊₁^3] * 
            [ϕαₖ₊₁ - ϕ₀ - dϕ₀*αₖ₊₁; ϕαₖ - ϕ₀ - dϕ₀*αₖ] / 
            (αₖ^2*αₖ₊₁^2*(αₖ₊₁-αₖ))
        αₖ₊₂ = (-b + sqrt(b^2 - 3a*dϕ₀)) / (3a)
        αₖ₊₂ = max(min(αₖ₊₂, αₖ₊₁ - ϵ), αₖ₊₁/2)
        αₖ = αₖ₊₁
        ϕαₖ = ϕαₖ₊₁
        αₖ₊₁ = T(αₖ₊₂)
    end
    return αₖ₊₁
end