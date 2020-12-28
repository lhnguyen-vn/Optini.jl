"""
    StaticInitial{T} <: AbstractInitial{T}

Static initial α₀ for line search methods. Newton and quasi-Newton algorithms for instance
should use α₀ = 1 as the initial guess.
"""
struct StaticInitial{T} <: AbstractInitial{T}
    α::T

    @doc """
    StaticInitial([α=1.0])

    Initialize `StaticInitial{T}`.
    """
    StaticInitial(α=1.0) = new{typeof(α)}(α)
end

(si::StaticInitial)(state, p) = fi.α

"""
    PreviousDecreaseInitial{T, D<:Ref} <: AbstractInitial{T}

`PreviousDecreaseInitial` assumes the first-order decrease is the same as the previous
iteration:

``α₀ = αₖ₋₁\\frac{∇f^\\intercal_{k-1} p_{k-1}}{∇f^\\intercal_k p_k}``

# Fields
- `α::T`: default step length for the first iteration
- `prev_decrease::D`: save the last iteration's first-order decrease
"""
struct PreviousDecreaseInitial{T, D<:Ref} <: AbstractInitial{T}
    α::T
    prev_decrease::D
end

"""
    PreviousDecreaseInitial([α::T=1.0]) where {T}

Initialize `PreviousDecreaseInitial`, with `prev_decrease` defaulting to `Ref(T(NaN))`.
"""
function PreviousDecreaseInitial(α::T=1.0) where {T}
    return PreviousDecreaseInitial(α, Ref(T(NaN)))
end

function reset!(pdi::PreviousDecreaseInitial)
    pdi.prev_decrease[] = NaN
    return nothing
end

function (pdi::PreviousDecreaseInitial{T})(state, p) where {T}
    prev_decrease = pdi.prev_decrease[]
    return isnan(prev_decrease) ? pdi.α : T(prev_decrease / (state.∇f ⋅ p))
end

function update!(pdi::PreviousDecreaseInitial, state, p, α)
    pdi.prev_decrease[] = α * (state.∇f ⋅ p)
    return nothing
end

"""
    QuadraticInitial{T} <: AbstractInitial{T}

`QuadraticInitial` assumes ``ϕ(α₀) - ϕ(0) = fₖ - fₖ₋₁`` and uses quadratic interpolation
with ``ϕ(α₀)``, ``ϕ(0)``, ``ϕ′(0)`` to compute the initial guess.

# Fields
- `α::T`: default step length for the first iteration
- `prev_f::Ref{T}`: save the last iteration's function value
"""
struct QuadraticInitial{T, D<:Ref} <: AbstractInitial{T}
    α::T
    prev_f::D
end

"""
    QuadraticInitial([α::T=1.0]) where {T}

Initialize `QuadraticInitial`, with `prev_f` defaulting to `Ref(T(NaN))`.
"""
function QuadraticInitial(α::T=1.0) where {T}
    return QuadraticInitial(α, Ref(T(NaN)))
end

function reset!(qi::QuadraticInitial)
    qi.prev_f[] = NaN
    return nothing
end

function (qi::QuadraticInitial{T})(state, p) where {T}
    prev_f = qi.prev_f[]
    return isnan(prev_f) ? qi.α : min(T(1.01*2*(state.f-qi.prev_f[])/(state.∇f⋅p)), one(T))
end

function update!(qi::QuadraticInitial, state, p, α)
    qi.prev_f[] = state.f
    return nothing
end