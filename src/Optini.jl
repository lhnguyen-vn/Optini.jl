module Optini

using LinearAlgebra
using Zygote

import Base.MathConstants: φ

# Generic types
include("types.jl")

# Bracketing algorithms
include("univariate/algorithms/fibonacci.jl")
include("univariate/algorithms/golden_section.jl")
include("univariate/algorithms/quadratic_fit.jl")
include("univariate/algorithms/bisection.jl")
include("univariate/bracket.jl")
include("univariate/optimize.jl")

# Line search methods
include("line_search/types.jl")
include("line_search/initial_step.jl")
include("line_search/algorithms/static.jl")
include("line_search/algorithms/exact.jl")
include("line_search/algorithms/backtracking.jl")
include("line_search/algorithms/interpolation.jl")
include("line_search/algorithms/strong_wolfe.jl")

# First order algorithms
include("multivariate/first_order/types.jl")
include("multivariate/first_order/algorithms/gradient_descent.jl")
include("multivariate/optimize.jl")

# Second order algorithms
include("multivariate/second_order/types.jl")
include("multivariate/second_order/algorithms/newton.jl")

export bracket, optimize
export Bisection, Fibonacci, GoldenSection, QuadraticFit

export StaticInitial, PreviousDecreaseInitial, QuadraticInitial
export StaticLineSearch, ExactLineSearch, BacktrackingLineSearch, InterpolationLineSearch,
    StrongWolfeLineSearch

export GradientDescent

export Newton

end