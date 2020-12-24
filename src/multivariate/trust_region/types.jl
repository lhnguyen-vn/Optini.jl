"""
    TrustRegion{H, M, D, R<:Ref, E} <: MultivariateAlgorithm

`TrustRegion` builds a quadratic model within a radius from the current point to approximate
the objective function and uses the model to choose the next step.

# Fields
- `hessian::H`: the method to compute or approximate the Hessian
- `method::M`: the method to compute the appropriate next step
- `Δ₀::D`: the initial radius
- `Δₘₐₓ::D`: the maximum radius
- `Δ::R`: a `Ref` containing the current radius
- `η::E`: the threshold to control when `TrustRegion` will take a step
"""
struct TrustRegion{H, M, D, E, S, R} <: MultivariateAlgorithm
    hessian::H
    method::M
    Δ₀::D
    Δₘₐₓ::D
    η::E
    ηₛ::E
    ηₑ::E
    σₛ::S
    σₑ::S
    Δ::R

    function TrustRegion(hessian, method, Δ₀::D, Δₘₐₓ::D, η, ηₛ, ηₑ, σₛ, σₑ) where {D}
        0 ≤ η < 0.5 || error("`η` must be in the interval [0, 0.5)")
        0 ≤ ηₛ < 0.5 || error("`ηₛ` must be in the interval [0, 0.5)")
        ηₑ ≥ 0.5 || error("`ηₑ` must be larger than or equal to 0.5")
        0 < σₛ < 1 || error("Shrinking factor `σₛ` must be in the interval (0, 1)")
        σₑ > 1 || error("Expanding factor `σₑ` must be larger than 1")
        Δ = Ref{D}()
        H = typeof(hessian)
        M = typeof(method)
        E = typeof(η)
        S = typeof(σₛ)
        R = typeof(Δ)
        new{H, M, D, E, S, R}(hessian, method, Δ₀, Δₘₐₓ, η, ηₛ, ηₑ, σₛ, σₑ, Δ)
    end
end

"""
    TrustRegion(; kwargs...)

Initiate `TrustRegion` algorithm.

# Keywords
- `hessian=Newton()`: the method to compute or approximate the Hessian
- `method=TwoDimSubspace()`: the method to compute the appropriate next step
- `delta=5.0`: the initial radius
- `deltamax=100.0`: the maximum radius
- `threshold=0.1`: the threshold to control when `TrustRegion` will take a step
- `shrinkthreshold=0.25`: the threshold to control when `TrustRegion` shrinks the current radius
- `expandthreshold=0.75`: the threshold to control when `TrustRegion` expands the current radius
- `shrinkfactor=0.25`: the factor used to shrink the trust region radius
- `expandfactor=2.0`: the factor used to expand the trust region radius
"""
function TrustRegion(; 
        hessian=Newton(), 
        method=TwoDimSubspace(), 
        delta=5.0, 
        deltamax=100.0, 
        threshold=0.1,
        shrinkthreshold=0.25,
        expandthreshold=0.75,
        shrinkfactor=0.25,
        expandfactor=2.0)
    D = promote_type(typeof(delta / 1), typeof(deltamax / 1))
    E = promote_type(typeof(threshold), typeof(shrinkthreshold), typeof(expandthreshold))
    S = promote_type(typeof(shrinkfactor), typeof(expandfactor))
    TrustRegion(hessian, method, D(delta), D(deltamax), E(threshold), E(shrinkthreshold), E(expandthreshold), 
        S(shrinkfactor), S(expandfactor))
end

order(tr::TrustRegion) = order(tr.hessian)