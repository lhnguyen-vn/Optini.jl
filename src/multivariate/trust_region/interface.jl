function reset!(alg::TrustRegion)
    alg.Δ[] = alg.Δ₀
    return nothing
end

function step!(x, alg::TrustRegion, f::Function, state)
    Δ = alg.Δ[]
    Δₘₐₓ = alg.Δₘₐₓ
    η = alg.η
    ηₛ = alg.ηₛ
    ηₑ = alg.ηₑ
    fx = state.f
    ∇fx = state.∇f
    B = approx_hessian(alg.hessian, state)
    model(p) = fx + ∇fx⋅p + (p'*B*p)/2
    p = alg.method(state, B, Δ)
    ρ = (fx - f(x + p)) / (model(zero(p)) - model(p))
    
    # Update radius
    if ρ ≥ ηₑ && norm(p) ≈ Δ # highly successful
        alg.Δ[] = min(alg.σₑ*Δ, Δₘₐₓ)
    elseif ρ < ηₛ
        alg.Δ[] = alg.σₛ*Δ
    end

    # Apply step
    if ρ ≥ η
        x .+= p
    end

    return nothing
end