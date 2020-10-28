"""
    InterpolationLineSearch{I<:AbstractInitial, C, T}

`InterpolationLineSearch` uses quadratic and cubic interpolations to compute the step length
satisfying the Armijo sufficient decrease condition.

# Fields
- `init::I`: the initial step length method
- `c::C`: the constant for the Armijo condition in the open interval (0, 1)
- `ϵ::T`: the minimum step length decrease from the initial step length
"""
struct InterpolationLineSearch{I<:AbstractInitial, C, T} <: AbstractLineSearch
    init::I
    c::C
    ϵ::T

    function InterpolationLineSearch(init, c, ϵ)
        0 < c < 1 || error("`c` must be in the open interval (0, 1)")
        I = typeof(init)
        C = typeof(c)
        T = typeof(ϵ)
        return new{I, C, T}(init, c, ϵ)
    end
end

"""
    InterpolationLineSearch(; kwargs...)

Initiate `InterpolationLineSearch`.

# Keywords
- `init=StaticInitial()`: the initial step length method
- `c=1e-4`: the constant for the Armijo condition in the open interval (0, 1)
- `ϵ=1e-6`: the minimum step length decrease from the initial step length
"""
function InterpolationLineSearch(; init=StaticInitial(), c=1e-4, ϵ=1e-6)
    InterpolationLineSearch(init, c, ϵ)
end

function (ils::InterpolationLineSearch{<:AbstractInitial{T}})(f, state, p) where {T}
    x = state.x
    ϕ₀ = state.f
    ∇y = state.∇f
    dϕ₀ = ∇y⋅p
    αₖ = ils.init(state, p)
    c = ils.c
    ϵ = ils.ϵ
    ϕ(α) = f(x + α * p)
    ϕαₖ = ϕ(αₖ)
    αₖ₊₁ = T(-(dϕ₀ * αₖ^2) / (2 * (ϕαₖ - ϕ₀ - dϕ₀ * αₖ)))
    while f(x + αₖ₊₁*p) > ϕ₀ + c*αₖ₊₁*dϕ₀
        ϕαₖ₊₁ = ϕ(αₖ₊₁)
        a, b = [αₖ^2 -αₖ₊₁^2; -αₖ^3 αₖ₊₁^3] * 
            [ϕαₖ₊₁ - ϕ₀ - dϕ₀*αₖ₊₁; ϕαₖ - ϕ₀ - dϕ₀*αₖ] / 
            (αₖ^2*αₖ₊₁^2*(αₖ₊₁-αₖ))
        α₂ = (-b + sqrt(b^2 - 3a*dϕ₀)) / (3a)
        α₂ = max(min(α₂, αₖ₊₁ - ϵ), αₖ₊₁/2)
        αₖ = αₖ₊₁
        αₖ₊₁ = T(α₂)
    end
    return αₖ₊₁
end