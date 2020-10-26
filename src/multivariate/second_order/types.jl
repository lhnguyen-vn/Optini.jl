struct SecondOrderState{Tx, Tf}
    x::Tx
    f::Tf
    ∇f::Vector{Tf}
    ∇²f::Matrix{Tf}
end

state(f, g, h, x, alg::SecondOrderAlgorithm) = SecondOrderState(x, f(x), g(x), h(x))