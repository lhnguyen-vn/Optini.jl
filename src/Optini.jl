module Optini

import MathOptInterface
const MOI = MathOptInterface

import Base.MathConstants: φ

# Generic method types
include("types.jl")

# Univariate optimizers
include("univariate/types.jl")
include("univariate/utils.jl")
include("univariate/optimize.jl")
include("univariate/algorithms/fibonacci.jl")
include("univariate/algorithms/golden_section.jl")
include("univariate/algorithms/quadratic_fit.jl")

export Fibonacci, GoldenSection, QuadraticFit
export optimize

end
