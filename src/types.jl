"""
    AbstractAlgorithm

Abstract supertype of all optimization algorithms.
"""
abstract type AbstractAlgorithm end

"""
    Solution{Tx, Tf}

Solution type for optimization models. 
"""
struct Solution{Tx, Tf}
    converged::Bool
    iter::Int
    minimizer::Tx
    minimum::Tf
end