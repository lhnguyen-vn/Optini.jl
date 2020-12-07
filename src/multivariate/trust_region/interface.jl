function reset!(alg::TrustRegion)
    alg.Δ[] = alg.Δ₀
    return nothing
end

function step!(x, alg::TrustRegion, f::Function, state)
    Δ = alg.Δ[]
    η = alg.η
    fx = state.f
    ∇fx = state.∇f
    B = approx_hessian(alg.hessian, state)
    model(p) = fx + ∇fx⋅p + (p'*B*p)/2
    p = alg.method(state, B, Δ)
    ρ = (fx - f(x + p)) / (model(zero(p)) - model(p))
    
    # Update radius
    if ρ < 0.25
        alg.Δ[] = Δ / 4
    else
        if ρ > 0.75 && norm(p) ≈ Δ
            alg.Δ[] = min(2Δ, alg.Δₘₐₓ)
        end
    end

    # Apply step
    if ρ > η
        x .+= p 
    end
end