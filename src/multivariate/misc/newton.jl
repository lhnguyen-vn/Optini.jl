"""
    Newton{M}

`Newton` dispatches to the Newton direction for `LineSearch` algorithm and uses the Hessian 
in `TrustRegion` models.

# Fields
- `modify::M`: the Hessian modification method to enforce positive definiteness
"""
struct Newton{M}
    modify::M
end

Newton() = Newton(MultipleIdentity())

order(::Newton) = SecondOrder()

function approx_hessian(n::Newton, state)
    ∇²f = state.∇²f
    isposdef(∇²f) && return ∇²f
    return n.modify(∇²f)
end

function direction(n::Newton, state)
    ∇f = state.∇f
    ∇²f = approx_hessian(n::Newton, state)
    return -inv(∇²f) * ∇f
end

"""
    MultipleIdentity{B, S}

Modify the Hessian by adding some multiples, τ, of the identity matrix, as described in 
Nocedal and Wright's Numerical Optimization.

# Fields
- `β::B`: constant to compute appropriate τ
- `scale::S`: scale factor to expand τ
"""
struct MultipleIdentity{B, S}
    β::B
    scale::S

    function MultipleIdentity(β, ρ)
        β > 0 || error("β must be positive")
        ρ > 1 || error("ρ must be larger than 1")
        new{typeof(β), typeof(ρ)}(β, ρ)
    end
end

"""
    MultipleIdentity(; beta=1e-3, scale=2.0)

Initiate `MultipleIdentity` method to modify the Hessian.
"""
MultipleIdentity(; beta=1e-3, scale=2.0) = MultipleIdentity(beta, scale)

function (mi::MultipleIdentity)(∇²f::Matrix{T}) where T
    β = T(mi.β)
    scale = T(mi.scale)
    τ = -minimum(diag(∇²f)) + β
    h = similar(∇²f)
    while true
        h .= ∇²f + τ*I
        isposdef(h) && return h
        τ = max(scale*τ, β)
    end
end