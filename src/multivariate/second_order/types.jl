struct SecondOrderState{Tx, Tf}
    x::Tx
    f::Tf
    ∇f::Vector{Tf}
    ∇²f::Matrix{Tf}
end

function SecondOrderState(f, x)
    fx, back = Zygote.pullback(f, x)
    ∇fx = back(one(typeof(fx)))[1]
    ∇²fx = Zygote.hessian(f, x)
    return SecondOrderState(x, fx, ∇fx, ∇²fx)
end

state(f, x, alg::SecondOrderAlgorithm) = SecondOrderState(f, x)