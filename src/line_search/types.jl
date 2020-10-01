function line_search(f, ∇f, x, d)
    objective(α) = f(x + α*d)
    solution = optimize(objective)
    if solution.converged
        return solution.minimizer
    end
end

function backtracking_line_search(f, ∇f, x, d, α=1.0; p=0.5, β=1e-4)
    y, g = f(x), ∇f(x)
    while f(x + α*d) > y + β*α*(g⋅d)
        α *= p
    end
    return α
end