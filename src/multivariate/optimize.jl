function optimize(f::Function, x::Vector{<:AbstractFloat};
        alg::MultivariateAlgorithm, abs_tol, max_iter, trace)
    s = state(f, x, alg)
    t = typeof(s)[]
    converged = false
    iter = 0
    while iter <= max_iter
        trace && push!(t, s)
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
        s = state(f, x, alg)
    end
    return Solution(converged, iter, s.x, s.f, Dict(:trace => t))
end

function update!(alg::MultivariateAlgorithm, state, p, α)
    update!(alg.linesearch, state, p, α)
    return alg
end