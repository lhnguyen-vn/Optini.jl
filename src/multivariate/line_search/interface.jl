function reset!(ai::AbstractInitial) end

function update!(ai::AbstractInitial, state, p, α) end

function reset!(alg::LineSearch)
    reset!(alg.initial)
end

function update!(alg, state, p, α)
    update!(alg.initial, state, p, α)
    return nothing
end

function step!(x, alg, f::Function, state)
    p = direction(alg.direction, state)
    α₀ = alg.initial(state, p)
    α = alg.method(f, state, p, α₀)
    update!(alg, state, p, α)
    x .+= α*p
end