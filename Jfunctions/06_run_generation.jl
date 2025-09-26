################# Function to run a single generation ###################

function runGeneration(;k, agents, mF, mC, nAgents, agentLength, testFunc, fCurrent, nMin, nMax, H = nothing, bounds = nothing, best = nothing, p = nothing, terminalValue = nothing, fEvals = nothing, maxFuncEvals = nothing)
    F = changeF(H = H, popSize = nAgents, mF = mF)
    C = changeC(H = H, popSize = nAgents, mC = mC, terminalValue = terminalValue)
    sF = Array{Float64}(undef, 0)
    sC = Array{Float64}(undef, 0)
    fDel = Array{Float64}(undef, 0)
    mut = deepcopy(agents)
    fDelMut = Array{Float64}(undef, nAgents)
    fMut = Array{Float64}(undef, nAgents)
    # for i in 1:nAgents
    #     mut[:,:,i], fDelMut[i], fMut[i] = mutateAgent(;iter = i, nAgents=nAgents, agents=agents, agentLength = agentLength, testFunc = testFunc, fCurrent = fCurrent[i], fVal = F, cVal = C, best = best, bounds = bounds)
    # end
    Threads.@threads for i in 1:nAgents
        local fi = fCurrent[i]
        local (m_slice, d, fm) = mutateAgent(;iter = i, nAgents=nAgents, agents=agents, agentLength = agentLength, testFunc = testFunc, fCurrent = fi, fVal = F, cVal = C, best = best, bounds = bounds)
        mut[:, :, i] = m_slice
        fDelMut[i] = d
        fMut[i] = fm
    end
    sC, sF, fDel, agents, fCurrent = updateParams(;mut = mut, fCurrent = fCurrent, fDelMut = fDelMut, fMut = fMut, nAgents=nAgents, agents=agents, fVal = F, cVal = C, sC = sC, sF = sF, fDel = fDel)
 
    mC[k] = updateMemC(fDel = fDel, mC = mC, sC = sC, k = k, terminalValue = terminalValue)
    mF[k] = updateMemF(fDel = fDel, mF = mF, sF = sF, k = k)

    k = updateK(k, H)
    keepAgents = cullHerd(nMin=nMin, nMax=nMax, maxFuncEvals = maxFuncEvals, fEvals = fEvals, fCurrent = fCurrent)
    agents = agents[:,:,keepAgents]
    nAgents = length(keepAgents)
    fCurrent = fCurrent[keepAgents]
    best = updateBest(p = p, nAgents = nAgents, fCurrent = fCurrent)
    return mC, mF, k, agents, nAgents, fCurrent, best
end
