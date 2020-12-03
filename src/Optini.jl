module Optini

using LinearAlgebra
using Zygote

import Base.MathConstants: φ

# Generic types
include("types.jl")

# Univariate algorithms
include("univariate/algorithms/fibonacci.jl")
include("univariate/algorithms/golden_section.jl")
include("univariate/algorithms/quadratic_fit.jl")
include("univariate/algorithms/bisection.jl")
include("univariate/bracket.jl")
include("univariate/optimize.jl")

# Multivariate algorithms
include("multivariate/types.jl")
include("multivariate/optimize.jl")
include("multivariate/misc/steepest.jl")
include("multivariate/misc/newton.jl")

# Line search methods
include("multivariate/line_search/types.jl")
include("multivariate/line_search/interface.jl")
include("multivariate/line_search/initial_guess.jl")
include("multivariate/line_search/algorithms/static.jl")
include("multivariate/line_search/algorithms/exact.jl")
include("multivariate/line_search/algorithms/backtracking.jl")
include("multivariate/line_search/algorithms/interpolation.jl")
include("multivariate/line_search/algorithms/strong_wolfe.jl")

# Trust region methods
include("multivariate/trust_region/types.jl")
include("multivariate/trust_region/interface.jl")
include("multivariate/trust_region/algorithms/cauchy.jl")
include("multivariate/trust_region/algorithms/dogleg.jl")
include("multivariate/trust_region/algorithms/2d_subspace.jl")

export bracket
export Bisection, Fibonacci, GoldenSection, QuadraticFit

export Steepest
export Newton

export StaticInitial, PreviousDecreaseInitial, QuadraticInitial
export StaticLineSearch, ExactLineSearch, BacktrackingLineSearch, InterpolationLineSearch,
    StrongWolfeLineSearch
export LineSearch, GradientDescent

export CauchyPoint, Dogleg, TwoDimSubspace
export TrustRegion

export optimize

end