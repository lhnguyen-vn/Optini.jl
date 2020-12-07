struct TwoDimSubspace end

function (tds::TwoDimSubspace)(state, B, Δ)
    isposdef(B) || 
        error("Two-dimensional subspace method only works for positive definite B.")
    p = similar(state.x)
    p1 = ∇f = state.∇f
    p2 = inv(B) * ∇f
    norm_∇f = norm(∇f)
    norm_p2 = norm(p2)
    prod = ∇f' * B * ∇f
    inv_prod = ∇f' * inv(B) * ∇f

    if (norm_p2 ≤ Δ) || (norm_p2 ≈ Δ) # analytical solution
        p .= -p2
    else if prod * inv_prod ≈ norm_∇f^4 # fall back to Cauchy Point
        _cauchypoint!(p, Δ, ∇f, norm_∇f, prod)
    else # solve constrained two-dimensional problem
        a = norm_∇f^2
        b = norm_p2^2
        c = prod
        d = inv_prod
        B̃ = [c a; a d]
        g̃ = [a; d]
        B̄ = [a d; d b]
        poly = λ -> (b - Δ^2)*(a^2 - c*d)^2 -
            2*(a^2 - c*d)*(Δ^2*a*d - Δ^2*b*c + a*b*d - d^3)*λ +
            (2*Δ^2*a^3*b - 3*Δ^2*a^2*d^2 - Δ^2*b^2*c^2 + 2*Δ^2*c*d^3 + a^3*b^2 - 2*a^2*b*d^2 + a*d^4)*λ^2 +
            2*Δ^2*(a*b - d^2)*(a*d - b*c)*λ^3 -
            Δ^2*(a*b - d^2)^2*λ^4
        λ = _positive_root(poly)
        α, β = -inv(B̃ + λ*B̄)*g̃
        @. p = α*p1 + β*p2
    end
    return p
end

function _positive_root(poly)
    step = 0.01
    scale = 2.0
    λ = step
    while poly(λ) > 0
        step *= scale
        λ += step
    end
    optimize(poly, 0.0, λ; alg=Bisection())
end