function optimize(f::Function, g::Function, h::Function, x::AbstractVector{T};
        alg::MultivariateAlgorithm, 
        abstol=1e-12, 
        maxiter::Integer=1_000_000, 
        savetrace::Bool=false,
        callback=()->nothing) where {T<:AbstractFloat}
    reset!(alg)
    x = copy(x)
    metadata = Dict{Symbol, Any}()
    curr_state = state(order(alg), f, g, h, x)
    trace = typeof(curr_state)[]
    converged = false
    iter = 0
    while iter <= maxiter
        callback()
        if savetrace
            push!(trace, curr_state)
            x = copy(x)
        end
        if norm(curr_state.∇f) < abstol
            converged = true
            break
        end
        (iter == maxiter || isinf(curr_state.f)) && break
        iter += 1
        step!(x, alg, f, curr_state)
        curr_state = state(order(alg), f, g, h, x)
    end
    savetrace && (metadata[:trace] = trace)
    return Solution(converged, iter, curr_state.x, curr_state.f, metadata)
end

function optimize(f::Function, g::Function, x; kwargs...)
    h = x -> Zygote.hessian(f, x)
    return optimize(f, g, h, x; kwargs...)
end

function optimize(f::Function, x; kwargs...)
    g = x -> Zygote.gradient(f, x)[1]
    h = x -> Zygote.hessian(f, x)
    return optimize(f, g, h, x; kwargs...)
end