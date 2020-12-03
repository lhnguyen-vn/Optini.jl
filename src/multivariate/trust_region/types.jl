"""
    TrustRegion{H, M, D, Dm, R<:Ref, E} <: MultivariateAlgorithm

`TrustRegion` builds a quadratic model within a radius from the current point to approximate
the objective function and uses the model to choose the next step.

# Fields
- `hessian::H`: the method to compute or approximate the Hessian
- `method::M<:AbstractTrustRegionMethod`: the method to compute the appropriate next step
- `Δ₀::D`: the initial radius
- `Δ::R`: a `Ref` containing the current radius
- `η::E`: the threshold to control when `TrustRegion` will take a step
"""
struct TrustRegion{H, M, D, Dm, R<:Ref, E} <: MultivariateAlgorithm
    hessian::H
    method::M
    Δ₀::D
    Δₘₐₓ::Dm
    η::E
    Δ::R

    function TrustRegion(hessian, method, Δ₀, Δₘₐₓ, η, Δ=Ref{typeof(Δ₀)}())
        0 ≤ η < 0.25 || error("`η` must be in the interval [0, 0.25)")
        H = typeof(hessian)
        M = typeof(method)
        D = typeof(Δ₀)
        Dm = typeof(Δₘₐₓ)
        E = typeof(η)
        R = typeof(Δ)
        TrustRegion{H, M, D, Dm, R, E}(hessian, method, Δ₀, Δₘₐₓ, η, Δ)
    end
end

"""
    TrustRegion(; kwargs...)

Initiate `TrustRegion` algorithm.

# Keywords
- `hessian`: the method to compute or approximate the Hessian
- `method`: the method to compute the appropriate next step
- `delta`: the initial radius
- `delta_max`: the maximum radius
- `eta`: the threshold to control when `TrustRegion` will take a step
"""
function TrustRegion(; hessian, method, delta, delta_max, eta)
    TrustRegion(hessian, method, delta, delta_max, eta)
end

order(tr::TrustRegion) = order(tr.hessian)