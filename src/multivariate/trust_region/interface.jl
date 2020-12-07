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
    model =  x -> fx + ∇fx⋅x + (x'*B*x)/2
    p = alg.method(state, B, Δ)
    new_x = x + p
    ρ = (fx - f(new_x)) / (model(x) - model(new_x))
    
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
        x .= new_x
    end
end