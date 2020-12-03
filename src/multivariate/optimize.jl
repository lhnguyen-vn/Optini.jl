function optimize(f::Function, g::Function, h::Function, x::Vector{T};
        alg::MultivariateAlgorithm, 
        abs_tol=1e-12, 
        max_iter::Integer=1_000_000, 
        save_trace::Bool=false,
        callback=()->nothing) where {T<:AbstractFloat}
    reset!(alg)
    metadata = Dict{Symbol, Any}()
    curr_state = state(order(alg), f, g, h, x)
    trace = typeof(curr_state)[]
    converged = false
    iter = 0
    while iter <= max_iter
        callback()
        if save_trace
            push!(trace, curr_state)
            x = copy(x)
        end
        if norm(curr_state.∇f) < abs_tol
            converged = true
            break
        end
        (iter == max_iter || isinf(curr_state.f)) && break
        iter += 1
        step!(x, alg, f, curr_state)
        curr_state = state(order(alg), f, g, h, x)
    end
    save_trace && (metadata[:trace] = trace)
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