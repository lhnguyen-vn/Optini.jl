abstract type AbstractInitial{T} end 

update!(ai::AbstractInitial, state, p, α) = ai

abstract type AbstractLineSearch end

function update!(ls::AbstractLineSearch, state, p, α)
    update!(ls.init, state, p, α)
    return ls
end