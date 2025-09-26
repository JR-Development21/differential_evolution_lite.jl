######################## Functions to update the population ###################

function sampPop(; iter, nAgents, best=nothing)
    pop_ind = iter
    pbest_indices = setdiff(best, iter)
    if length(pbest_indices) > 1
        pbest_ind = sample(pbest_indices, 1, replace = false)
    else
        pbest_ind = pbest_indices
    end
    r1_ind = sample(setdiff(1:nAgents, iter), 1, replace = false)
    r2_ind = sample(setdiff(1:nAgents, [iter, r1_ind]), 1, replace = false)

    pop_ind = vcat(pop_ind, pbest_ind, r1_ind, r2_ind)
    return pop_ind
end

function updateBest(; p = nothing, nAgents = nothing, fCurrent = nothing)

    order = sortperm(fCurrent)
    p_num = round(maximum([2, p * nAgents]))
    p_num = Int(p_num)
    pbest_ind = Array{Float64}(undef, p_num)
    pbest_ind = order[1:p_num]
    return pbest_ind
end

function cullHerd(; nMin = nothing, nMax = nothing, maxFuncEvals = nothing, fEvals = nothing, fCurrent = nothing)
    nAgents = Int(round((nMin - nMax) / maxFuncEvals * fEvals + nMax))
    keepAgents = sortperm(fCurrent)[1:nAgents]
    return keepAgents
end