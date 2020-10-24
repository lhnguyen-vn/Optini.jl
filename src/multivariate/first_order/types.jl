struct FirstOrderState{Tx, Tf}
    x::Tx
    f::Tf
    ∇f::Vector{Tf}
end

function FirstOrderState(f, x)
    fx, back = Zygote.pullback(f, x)
    ∇fx = back(one(typeof(fx)))[1]
    return FirstOrderState(x, fx, ∇fx)
end

state(f, x, alg::FirstOrderAlgorithm) = FirstOrderState(f, x)