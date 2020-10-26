abstract type AbstractInitial{T} end

reset!(ai::AbstractInitial) = ai

update!(ai::AbstractInitial, state, p, α) = ai

abstract type AbstractLineSearch end

function reset!(ls::AbstractLineSearch)
    reset!(ls.init)
    return ls
end

function update!(ls::AbstractLineSearch, state, p, α)
    update!(ls.init, state, p, α)
    return ls
end