"""
    AbstractAlgorithm

Abstract supertype of all optimization algorithms.
"""
abstract type AbstractAlgorithm end

"""
    UnivariateAlgorithm

Abstract supertype of univariate optimization algorithms.
"""
abstract type UnivariateAlgorithm <: AbstractAlgorithm end

"""
    MultivariateAlgorithm

Abstract supertype of multivariate optimization algorithms.
"""
abstract type MultivariateAlgorithm <: AbstractAlgorithm end

"""
    Solution{Tx, Tf}

Solution type for optimization models. 
"""
struct Solution{Tx, Tf, M}
    converged::Bool
    iter::Int
    minimizer::Tx
    minimum::Tf
    metadata::M
end

Solution(converged, iter, minimizer, minimum) = Solution(converged, iter, minimizer, minimum, Dict())