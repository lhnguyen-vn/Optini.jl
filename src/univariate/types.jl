struct UnivariateSolution{Tx, Tf} <: AbstractSolution
    converged::Bool
    iter::Int
    minimizer::Tx
    minimum::Tf
end