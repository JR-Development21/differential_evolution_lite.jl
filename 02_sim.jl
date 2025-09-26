t1 = time()
cd("C:\\Users\\A02180741\\Documents\\differential_evolution.jl\\")

using Distributed
length(Sys.cpu_info())    # Number of logical CPU cores available in the system
nprocs()      # Number of available processes.
nworkers()

#
addprocs(min(length(Sys.cpu_info()) - 3, 30))
nprocs()
nworkers()

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
    nReps = 140

    ## set up  combinatoric runs by factor numAgents
    nVec = vcat(3:9, 6:12, 10:16)
    kVec = repeat(1:3, inner = 7)
    agentVec = [50, 150, 500]
    nRuns = length(kVec)
    nSize = length(agentVec)
    critVec = [D_criterion, A_criterion, G_criterion, I_criterion]
    totalRuns = nRuns * nReps * nSize
    DVec = Array{Float64}(undef, totalRuns, 6)
    AVec = Array{Float64}(undef, totalRuns, 6)
    GVec = Array{Float64}(undef, totalRuns, 6)
    IVec = Array{Float64}(undef, totalRuns, 6)

    DDesign = Array{Array{Float64}}(undef, totalRuns)
    ADesign = Array{Array{Float64}}(undef, totalRuns)
    GDesign = Array{Array{Float64}}(undef, totalRuns)
    IDesign = Array{Array{Float64}}(undef, totalRuns)


    function run_optimization(i, j, k, l)
        println("Combo = ", i, ", nrep = ", j)
        K = kVec[i]
        N = nVec[i]
        d = K * N
        bounds = easyBounds(K)
        nAgents = agentVec[l]
        agents = initializeAgents(nAgents, N, K, bounds)
        pArchive = 1
        nArchive = Int(floor(nAgents * pArchive))
        archive = initializeAgents(nArchive, N, K, bounds)
        bounds = [[-1, 1] for i in 1:d]
        C = [.8 for x in 1:nAgents]
        F = [.7 for x in 1:nAgents]
        H = 6
        mF = [.5 for x in 1:H]
        mC = [.5 for x in 1:H]
        crit = x -> critVec[k](x; N = N, K = K)
        p = .11

        out = differentialEvolution(deepcopy(agents), C, F, crit, maxItersPSO = 5000, convergeIters = 100, mF = mF, mC = mC, method = "current_to_pbest/1", variant = "SHADE", p = p, H = H, doArchive = true, archive = archive, nArchive = nArchive, pArchive = pArchive, bounds = bounds, boundType = "absorbingWall")

        if k == 1
            critVal = D_criterion(out[1]; N = N, K = K)
            return (:D, j + nReps * nSize * (i - 1) + nReps * (l - 1), [j, K, N, out[3], nAgents, critVal], out[1])
        elseif k == 2
            critVal = A_criterion(out[1]; N = N, K = K)
            return (:A, j + nReps * nSize * (i - 1) + nReps * (l - 1), [j, K, N, out[3], nAgents, critVal], out[1])
        elseif k == 3
            critVal = G_criterion(out[1]; N = N, K = K)
            return (:G, j + nReps * nSize * (i - 1) + nReps * (l - 1), [j, K, N, out[3], nAgents, critVal], out[1])
        else
            critVal = I_criterion(out[1]; N = N, K = K)
            return (:I, j + nReps * nSize * (i - 1) + nReps * (l - 1), [j, K, N, out[3], nAgents, critVal], out[1])
        end
    end
end

## Use pmap to parallelize the optimization runs
results = pmap(x -> run_optimization(x[1], x[2], x[3], x[4]), [(i, j, k, l) for i in 1:nRuns for j in 1:nReps for k in 1:4 for l in 1:nSize])

## Store the results in the respective matrices
for res in results
    critType, index, vec, design = res
    if critType == :D
        DVec[index, :] = vec
        DDesign[index] = design
    elseif critType == :A
        AVec[index, :] = vec
        ADesign[index] = design
    elseif critType == :G
        GVec[index, :] = vec
        GDesign[index] = design
    else
        IVec[index, :] = vec
        IDesign[index] = design
    end
end


D_frame = DataFrame(hcat(DVec, DDesign), ["testRun","k","N","niter","Sv","fgbest","gbest"])
CSV.write("D_sim.csv", D_frame)

A_frame = DataFrame(hcat(AVec, ADesign), ["testRun","k","N","niter","Sv","fgbest","gbest"])
CSV.write("A_sim.csv", A_frame)

G_frame = DataFrame(hcat(GVec, GDesign), ["testRun","k","N","niter","Sv","fgbest","gbest"])
CSV.write("G_sim.csv", G_frame)

I_frame = DataFrame(hcat(IVec, IDesign), ["testRun","k","N","niter","Sv","fgbest","gbest"])
CSV.write("I_sim.csv", I_frame)

elapsed_time = time() - t1
println("Elapsed time: ", elapsed_time, " seconds")