struct GradientDescent{T} <: FirstOrderAlgorithm
    line_search::T
end

GradientDescent(α::T) where {T<:AbstractFloat} = GradientDescent((f, ∇f, x, d) -> α)

function optimize(f::Function, x, algorithm::GradientDescent; 
        rel_tol, 
        abs_tol, 
        max_iter, 
        in_place=false)
    if !in_place
        x = copy(x)
    end
    converged = false
    iter = 1
    while iter < max_iter
        ∇f(x) = Zygote.gradient(f, x)[1]
        d = -∇f(x)
        α = algorithm.line_search(f, ∇f, x, d)
        α === nothing && throw(ErrorException("No step length found"))
        step = α * d
        x_tol = x * rel_tol .+ abs_tol
        if norm(step) < norm(x_tol)
            converged = true
            break
        end
        x .+= step
        iter += 1
    end
    return Solution(converged, iter, x, f(x))
end