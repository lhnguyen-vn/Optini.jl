"""
    Steepest

`Steepest` dispatches to the opposite of the gradient direction for `LineSearch`.
"""
struct Steepest end

order(::Steepest) = FirstOrder()

direction(::Steepest, state) = -state.∇f