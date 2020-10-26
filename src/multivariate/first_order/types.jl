struct FirstOrderState{Tx, Tf}
    x::Tx
    f::Tf
    ∇f::Vector{Tf}
end

state(f, g, h, x, alg::FirstOrderAlgorithm) = FirstOrderState(x, f(x), g(x))