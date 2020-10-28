"""
    ExactLineSearch{I<:AbstractInitial, A<:UnivariateAlgorithm, O1, O2} <: AbstractLineSearch

`ExactLineSearch` computes the optimal step length using univariate methods. Since most 
univariate algorithms require repeated calculations of the objective to shrink the bracket
around a local minimum, `ExactLineSearch` is expected to be expensive.

# Fields
- `init::I`: the initial step length method
- `alg::A`: the univariate algorithm
- `bracket_options::O1`: named tuple of bracketing options
- `alg_options::O2`: named tuple of the univariate algorithm options
"""
struct ExactLineSearch{I<:AbstractInitial, A<:UnivariateAlgorithm, O1, O2} <: AbstractLineSearch
    init::I
    alg::A
    bracket_options::O1
    alg_options::O2
end

"""
    ExactLineSearch(; kwargs...)

Initiate a `ExactLineSearch`. Use keyword arguments to tweak the bracketing phase and the 
univariate algorithm.

# Keywords
- `init=StaticInitial()`: the initial step length method
- `step=0.01`: bracketing step size
- `scale=2.0`: bracketing scale factor
- `bracket_max_iter::Integer=100`: bracketing maximum number of iterations
- `alg::UnivariateAlgorithm`: the univariate algorithm to optimize the step length
- `rel_tol=1e-6`: relative tolerance
- `abs_tol=1e-12`: absolute tolerance
- `max_iter::Integer=1_000`: maximum number of iterations
"""
function ExactLineSearch(;
        init=StaticInitial(),
        step=0.01,
        scale=2.0,
        bracket_max_iter::Integer=100,
        alg::UnivariateAlgorithm=GoldenSection(),
        rel_tol=1e-6,
        abs_tol=1e-12,
        max_iter::Integer=1_000)
    alg_options = (; rel_tol, abs_tol, max_iter)
    bracket_options = (; step, scale, max_iter=bracket_max_iter)
    ExactLineSearch(init, alg, bracket_options, alg_options)
end

reset!(els::ExactLineSearch) = els

function (els::ExactLineSearch)(f, state, p)
    α₀ = els.init(state, p)
    ϕ(α) = f(state.x + α * p)
    br = bracket(ϕ, α₀; els.bracket_options...)
    isnothing(br) && error("Initial minimum bracketing failed.")
    lower, upper = br
    solution = optimize(ϕ, lower, upper; els.alg, els.alg_options...)
    solution.converged || error("Exact line search failed.")
    return solution.minimizer
end

update!(els::ExactLineSearch, state, p, α) = els