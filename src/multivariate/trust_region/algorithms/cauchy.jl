struct CauchyPoint end

function (cp::CauchyPoint)(state, B, Δ)
    p = similar(state.x)
    ∇f = state.∇f
    norm_∇f = norm(∇f)
    p .= -Δ/norm_∇f*∇f
    prod = ∇f'*B*∇f
    if prod > 0
        p .*= min(1, norm_∇f^3/(Δ*prod))
    end
    return p
end