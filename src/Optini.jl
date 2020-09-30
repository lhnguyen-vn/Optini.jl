module Optini

using Zygote

import MathOptInterface
const MOI = MathOptInterface

import Base.MathConstants: φ

# Generic method types
include("types.jl")

# bracketing optimizers
include("bracketing/types.jl")
include("bracketing/utils.jl")
include("bracketing/optimize.jl")
include("bracketing/algorithms/fibonacci.jl")
include("bracketing/algorithms/golden_section.jl")
include("bracketing/algorithms/quadratic_fit.jl")
include("bracketing/algorithms/bisection.jl")

export Fibonacci, GoldenSection, QuadraticFit, Bisection
export optimize

end
