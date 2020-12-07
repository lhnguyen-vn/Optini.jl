abstract type AbstractOrder end

struct FirstOrder <: AbstractOrder end

struct SecondOrder <: AbstractOrder end

struct FirstOrderState{Tx, Tf}
    x::Tx
    f::Tf
    ∇f::Vector{Tf}
end

state(::FirstOrder, f, g, h, x) = FirstOrderState(x, f(x), g(x))

struct SecondOrderState{Tx, Tf}
    x::Tx
    f::Tf
    ∇f::Vector{Tf}
    ∇²f::Matrix{Tf}
end

state(::SecondOrder, f, g, h, x) = SecondOrderState(x, f(x), g(x), h(x))