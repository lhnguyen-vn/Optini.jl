struct Dogleg end

function (dl::Dogleg)(state, B, Δ)
    p = similar(state.x)
    ∇f = state.∇f
    norm_∇f = norm(∇f)
    prod = ∇f' * B * ∇f
    pᶜ = -norm_∇f^2 / prod*∇f
    norm_pᶜ = norm(pᶜ)
    if Δ ≤ norm_pᶜ
        p .= Δ / norm_pᶜ * pᶜ
    else
        isposdef(B) || error("Dogleg method only works for positive definite B")
        pᴮ = -inv(B)*∇f
        if Δ ≥ pᴮ
            p .= pᴮ
        else
            pᴮpᶜ = pᴮ⋅pᶜ
            norm_pᴮ = norm(pᴮ)
            r = (norm_pᴮ^2 - pᴮpᶜ) / (norm_pᶜ^2 + norm_pᴮ^2 - 2pᴮpᶜ)
            @. p = r*pᶜ + (1-r)*pᴮ
        end
    end
    return p
end