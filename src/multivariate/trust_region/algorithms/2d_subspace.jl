struct TwoDimSubspace end

function (tds::TwoDimSubspace)(state, B, Δ)
    isposdef(B) || 
        error("Two-dimensional subspace method only works for positive definite B.")
    p = similar(state.x)
    p1 = ∇f = state.∇f
    p2 = inv(B) * ∇f
    norm_∇f = norm(∇f)
    prod = ∇f' * B * ∇f
    inv_prod = ∇f' * inv(B) * ∇f
    g̃ = [norm_∇f^2; inv_prod]
    B̃ = [prod norm_∇f^2; norm_∇f^2 inv_prod]
    α, β = -inv(B̃)*g̃
    @. p = α * p1 + β * p2
    if norm(p) > Δ^2
        
    end
    return p
end