abstract type BracketingAlgorithm <: AbstractAlgorithm end

struct BracketingSolution{Tx, Tf} <: AbstractSolution
    converged::Bool
    iter::Int
    minimizer::Tx
    minimum::Tf
end