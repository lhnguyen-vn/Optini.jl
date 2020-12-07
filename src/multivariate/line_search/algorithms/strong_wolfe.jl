"""
    StrongWolfeLineSearch{C1, C2, S, T}

`StrongWolfeLineSearch` computes a step length that satifies both the Armijo sufficient 
decrease and the curvature conditions.

# Fields
- `c₁::C1`: Armijo condition constant
- `c₂::C2`: curvature condition constant
- `scale::S`: scale factor to expand step length bracket
- `α_max::T`: maximum step length
"""
struct StrongWolfeLineSearch{C1, C2, S, T}
    c₁::C1
    c₂::C2
    scale::S
    αₘₐₓ::T

    function StrongWolfeLineSearch(c₁, c₂, scale, αₘₐₓ)
        0 < c₁ < c₂ < 1 || error("Constants must satisfy 0 < `c₁` < `c₂` < 1")
        scale > 1 || error("`scale` factor must be larger than 1")
        C1 = typeof(c₁)
        C2 = typeof(c₂)
        S = typeof(scale)
        T = typeof(αₘₐₓ)
        return new{C1, C2, S, T}(c₁, c₂, scale, αₘₐₓ)
    end
end

"""
    StrongWolfeLineSearch(; kwargs...)

Initiate `StrongWolfeLineSearch`.

# Keywords
- `c1=1e-4`: Armijo condition constant
- `c2=0.9: curvature condition constant
- `scale=2.0`: scale factor to expand step length bracket
- `alphamax=100_000`: maximum step length
"""
function StrongWolfeLineSearch(;
        c1=1e-4, 
        c2=0.9,
        scale=2.0,
        alphamax=100_000)
    return StrongWolfeLineSearch(c1, c2, scale, alphamax)
end

function (sw::StrongWolfeLineSearch)(f, state, p, α₀::T) where {T}
    x = state.x
    ϕ₀ = state.f
    dϕ₀ = state.∇f ⋅ p
    c₁ = sw.c₁
    c₂ = sw.c₂
    scale = T(sw.scale)
    αₘₐₓ = T(sw.αₘₐₓ)
    αₖ = zero(T)
    ϕαₖ = ϕ₀
    αₖ₊₁ = α₀
    ϕ(α) = f(x + α*p)
    iter = 0
    while αₖ₊₁ < αₘₐₓ
        ϕαₖ₊₁ = ϕ(αₖ₊₁)
        if ϕαₖ₊₁ > ϕ₀ + c₁*αₖ₊₁*dϕ₀ || (ϕαₖ₊₁ ≥ ϕαₖ && iter > 1)
            return _zoom(ϕ, αₖ, αₖ₊₁, ϕ₀, dϕ₀, ϕαₖ, c₁, c₂)
        end
        dϕαₖ₊₁ = Zygote.gradient(ϕ, αₖ₊₁)[1]
        abs(dϕαₖ₊₁) ≤ -c₂*dϕ₀ && return αₖ₊₁
        if dϕαₖ₊₁ ≥ 0
            return _zoom(ϕ, αₖ₊₁, αₖ, ϕ₀, dϕ₀, ϕαₖ₊₁, c₁, c₂)
        end
        αₖ = αₖ₊₁
        ϕαₖ = ϕαₖ₊₁
        αₖ₊₁ *= scale
        iter += 1
    end
    return αₘₐₓ
end

function _zoom(ϕ, αₗₒ, αₕᵢ, ϕ₀, dϕ₀, ϕαₗₒ, c₁, c₂)
    while true
        α = (αₗₒ + αₕᵢ) / 2
        ϕα = ϕ(α)
        if ϕα > ϕ₀ + c₁*α*dϕ₀ || ϕα ≥ ϕαₗₒ
            αₕᵢ = α
        else
            dϕα = Zygote.gradient(ϕ, α)[1]
            abs(dϕα) ≤ -c₂*dϕ₀ && return α
            if dϕα*(αₕᵢ - αₗₒ) ≥ 0
                αₕᵢ = αₗₒ
            end
            αₗₒ = α
            ϕαₗₒ = ϕα
        end 
    end
end