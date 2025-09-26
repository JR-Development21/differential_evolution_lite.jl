############## Function to mutate a given agent ##############

function mutateAgent(;iter, nAgents, agents, agentLength, fVal, cVal, testFunc, fCurrent, best = nothing, bounds = nothing)
    popInd = sampPop(;iter=iter, nAgents = nAgents, best = best)

    iterChange = rand(1:agentLength)
    iterProb = rand(Uniform(), agentLength)
    iterProb[iterChange] = 0
    mutVec = deepcopy(agents[:,:,popInd[1]])
    for j in 1:agentLength
        if iterProb[j] <= cVal[1]
            mutVec[j] += fVal[1] * (agents[:,:,popInd[2]][j] - agents[:,:,popInd[1]][j]) + fVal[1] * (agents[:,:,popInd[3]][j] - agents[:,:,popInd[4]][j])
        end
    end


    for j in 1:agentLength
        if mutVec[j] < bounds[j][1]
            mutVec[j] = bounds[j][1]
        elseif mutVec[j] > bounds[j][2]
            mutVec[j] = bounds[j][2]
        end
    end

    fMut = testFunc(mutVec)
    fDelMut = fMut - fCurrent
    return mutVec, fDelMut, fMut
end

function updateParams(;agents, fCurrent, mut, fMut, fDelMut, sC, sF, fDel, cVal, fVal, nAgents)
    for i in 1:nAgents
        if fDelMut[i] < 0
            if ~isinf(fDelMut[i])
                sC = vcat(sC, cVal[i])
                sF = vcat(sF, fVal[i])
                fDel = vcat(fDel, abs(fDelMut[i]))
            end
            agents[:,:,i] = mut[:,:,i]
            fCurrent[i] = fMut[i]
        end
    end
    return sC, sF, fDel, agents, fCurrent
end
