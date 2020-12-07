"""
    ExactLineSearch{A<:UnivariateAlgorithm, O1, O2}

`ExactLineSearch` computes the optimal step length using univariate methods. Since most 
univariate algorithms require repeated calculations of the objective to shrink the bracket
around a local minimum, `ExactLineSearch` is expected to be expensive.

# Fields
- `alg::A`: the univariate algorithm
- `bracket_options::O1`: named tuple of bracketing options
- `alg_options::O2`: named tuple of the univariate algorithm options
"""
struct ExactLineSearch{A<:UnivariateAlgorithm, O1, O2}
    alg::A
    bracket_options::O1
    alg_options::O2
end

"""
    ExactLineSearch(; kwargs...)

Initiate a `ExactLineSearch`. Use keyword arguments to tweak the bracketing phase and the 
univariate algorithm.

# Keywords
- `alg::UnivariateAlgorithm`: the univariate algorithm to optimize the step length
- `step=0.01`: bracketing step size
- `scale=2.0`: bracketing scale factor
- `bracket_maxiter::Integer=100`: bracketing maximum number of iterations
- `reltol=1e-6`: relative tolerance
- `abstol=1e-12`: absolute tolerance
- `maxiter::Integer=1_000`: maximum number of iterations
"""
function ExactLineSearch(;
        alg::UnivariateAlgorithm=GoldenSection(),
        step=0.01,
        scale=2.0,
        bracket_maxiter::Integer=100,
        reltol=1e-6,
        abstol=1e-12,
        maxiter::Integer=1_000)
    alg_options = (; reltol, abstol, maxiter)
    bracket_options = (; step, scale, maxiter=bracket_maxiter)
    ExactLineSearch(alg, bracket_options, alg_options)
end

function (els::ExactLineSearch)(f, state, p, α₀)
    ϕ(α) = f(state.x + α * p)
    br = bracket(ϕ, α₀; els.bracket_options...)
    isnothing(br) && error("Initial minimum bracketing failed.")
    lower, upper = br
    solution = optimize(ϕ, lower, upper; els.alg, els.alg_options...)
    solution.converged || error("Exact line search failed.")
    return solution.minimizer
end