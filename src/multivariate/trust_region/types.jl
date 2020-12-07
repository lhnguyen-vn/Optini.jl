"""
    TrustRegion{H, M, D, R<:Ref, E} <: MultivariateAlgorithm

`TrustRegion` builds a quadratic model within a radius from the current point to approximate
the objective function and uses the model to choose the next step.

# Fields
- `hessian::H`: the method to compute or approximate the Hessian
- `method::M<:AbstractTrustRegionMethod`: the method to compute the appropriate next step
- `Δ₀::D`: the initial radius
- `Δ::R`: a `Ref` containing the current radius
- `η::E`: the threshold to control when `TrustRegion` will take a step
"""
struct TrustRegion{H, M, D, R<:Ref, E} <: MultivariateAlgorithm
    hessian::H
    method::M
    Δ₀::D
    Δₘₐₓ::D
    η::E
    Δ::R

    function TrustRegion(hessian, method, Δ₀::D, Δₘₐₓ::D, η, Δ=Ref{D}()) where D
        0 ≤ η < 0.25 || error("`η` must be in the interval [0, 0.25)")
        H = typeof(hessian)
        M = typeof(method)
        E = typeof(η)
        new{H, M, D, Ref{D}, E}(hessian, method, Δ₀, Δₘₐₓ, η, Δ)
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
- `eta=0.1`: the threshold to control when `TrustRegion` will take a step
"""
function TrustRegion(; 
        hessian=Newton(), 
        method=TwoDimSubspace(), 
        delta=5.0, 
        deltamax=100.0, 
        eta=0.1)
    T = promote_type(typeof(delta / 1), typeof(deltamax / 1))
    TrustRegion(hessian, method, T(delta), T(deltamax), eta)
end

order(tr::TrustRegion) = order(tr.hessian)