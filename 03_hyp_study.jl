cd("C:\\Users\\A02180741\\Documents\\differential_evolution.jl\\")

using Distributed, Dates
length(Sys.cpu_info())    # Number of logical CPU cores available in the system
nprocs()      # Number of available processes.
nworkers()

#
addprocs(min(length(Sys.cpu_info()) - 1, 30))
nprocs()
nworkers()

print(Dates.format(now(), "HH:MM"))

@everywhere using Distributions, LinearAlgebra, Combinatorics, Random, CSV, DataFrames, Distributed

@everywhere begin
    functions = readdir("Jfunctions/")
    functions = functions[findall(occursin.(".jl", functions))]
    nf = length(functions)
    for i in 1:nf
        ftemp = string("Jfunctions\\", functions[i])
        include(ftemp)
    end


    # Random.seed!(16)
    nReps = 100

    # create design from Walsh results
    hpoDesign = [1.0 -0.9999999999999996; 7.2761183106552234e-9 -1.5297205054337253e-8; -3.142222040344605e-8 -1.3101472581070075e-8; -1.0 -1.0; -6.145091794900058e-9 -0.9999999999999999; -1.6246069449272817e-8 1.0; 2.2726978421678407e-9 -3.953940820479788e-9; -4.718212683873243e-10 3.169048643833188e-9; -1.0 1.0; 1.0 -5.516752889729785e-9; 0.9999999999999996 1.0; -1.0 6.610875494983516e-10]
    rangePop = [50, 500]
    rangeP = [0, .2]
    agentVec = Int.(round.(design2Space(hpoDesign[:,1], rangePop[1], rangePop[2])))
    pVec = design2Space(hpoDesign[:,1], rangeP[1], rangeP[2])

    ## set up  combinatoric runs by factor numAgents
    nVec = [6, 12, 10, 16]
    kVec = [2, 2, 3, 3]
    nRuns = length(kVec)
    nSize = length(agentVec)
    critVec = [A_criterion, I_criterion]
    nCrit = length(critVec)
    totalRuns = nRuns * nReps * nSize
    AVec = Array{Float64}(undef, totalRuns, 6)
    IVec = Array{Float64}(undef, totalRuns, 6)

    ADesign = Array{Array{Float64}}(undef, totalRuns)
    IDesign = Array{Array{Float64}}(undef, totalRuns)


    function run_optimization(i, j, k, l)
        K = kVec[i]
        N = nVec[i]
        d = K * N
        bounds = easyBounds(K)
        nAgents = agentVec[l]
        agents = initializeAgents(nAgents, N, K, bounds)
        pArchive = 2.6
        nArchive = Int(floor(nAgents * pArchive))
        archive = initializeAgents(nArchive, N, K, bounds)
        bounds = [[-1, 1] for i in 1:d]
        C = [.8 for x in 1:nAgents]
        F = [.7 for x in 1:nAgents]
        H = 6
        mF = [.5 for x in 1:H]
        mC = [.5 for x in 1:H]
        crit = x -> critVec[k](x; N = N, K = K)
        p = pVec[l]

        out = differentialEvolution(deepcopy(agents), C, F, crit, maxItersPSO = 5000, convergeIters = 100, mF = mF, mC = mC, method = "current_to_pbest/1", variant = "L-SHADE", p = p, H = H, doArchive = true, archive = archive, nArchive = nArchive, pArchive = pArchive, bounds = bounds, boundType = "absorbingWall")


        if k == 1
            critVal = A_criterion(out[1]; N = N, K = K)
            return (:A, j + nReps * nSize * (i - 1) + nReps * (l - 1), [j, K, N, out[3], l, critVal], out[1])
        else
            critVal = I_criterion(out[1]; N = N, K = K)
            return (:I, j + nReps * nSize * (i - 1) + nReps * (l - 1), [j, K, N, out[3], l, critVal], out[1])
        end
    end
end

## Use pmap to parallelize the optimization runs
results = pmap(x -> run_optimization(x[1], x[2], x[3], x[4]), [(i, j, k, l) for i in 1:nRuns for j in 1:nReps for k in 1:nCrit for l in 1:nSize])

## Store the results in the respective matrices
for res in results
    critType, index, vec, design = res
    if critType == :A
        AVec[index, :] = vec
        ADesign[index] = design
    else
        IVec[index, :] = vec
        IDesign[index] = design
    end
end



A_frame = DataFrame(hcat(AVec, ADesign), ["testRun","k","N","niter","designRow","fgbest","gbest"])
CSV.write("A_sim_HPO.csv", A_frame)


I_frame = DataFrame(hcat(IVec, IDesign), ["testRun","k","N","niter","designRow","fgbest","gbest"])
CSV.write("I_sim_HPO.csv", I_frame)

print(Dates.format(now(), "HH:MM"))
