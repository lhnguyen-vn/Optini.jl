function optimize(f::Function, g::Function, h::Function, x::Vector{T};
        alg::MultivariateAlgorithm, 
        abs_tol=1e-12, 
        max_iter::Integer=1_000_000, 
        save_trace::Bool=false,
        callback=()->nothing) where {T<:AbstractFloat}
    reset!(alg)
    metadata = Dict()
    s = state(f, g, h, x, alg)
    t = typeof(s)[]
    converged = false
    iter = 0
    while iter <= max_iter
        callback()
        save_trace && push!(t, s)
        if norm(s.∇f) < abs_tol
            converged = true
            break
        end
        (iter == max_iter || isinf(s.f)) && break
        iter += 1
        p = alg(s)
        α = alg.linesearch(f, s, p)
        update!(alg, s, p, α)
        step = convert(typeof(x), α*p)
        x += step
        s = state(f, g, h, x, alg)
    end
    save_trace && (metadata[:trace] = t)
    return Solution(converged, iter, s.x, s.f, metadata)
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

function reset!(alg::MultivariateAlgorithm)
    reset!(alg.linesearch)
    return alg
end

function update!(alg::MultivariateAlgorithm, state, p, α)
    update!(alg.linesearch, state, p, α)
    return alg
end