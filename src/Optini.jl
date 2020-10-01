module Optini

using Zygote

import MathOptInterface
const MOI = MathOptInterface

import Base.MathConstants: φ

# Generic types
include("types.jl")

# Bracketing algorithms
include("bracketing/types.jl")
include("bracketing/utils.jl")
include("bracketing/optimize.jl")
include("bracketing/algorithms/fibonacci.jl")
include("bracketing/algorithms/golden_section.jl")
include("bracketing/algorithms/quadratic_fit.jl")
include("bracketing/algorithms/bisection.jl")

export Bisection, Fibonacci, GoldenSection, QuadraticFit
export bracket, optimize

end
