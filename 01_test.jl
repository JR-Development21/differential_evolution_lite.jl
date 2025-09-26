using Distributions, LinearAlgebra, Distances, Plots, Base.Threads

cd("C:/Users/ritch/Documents/DOE Papers/differential_evolution_lite.jl/")

functions = readdir("Jfunctions/")
functions = functions[findall(occursin.(".jl", functions))]
nf = length(functions)
for i in 1:nf
    ftemp = string("Jfunctions\\", functions[i])
    include(ftemp)
end

using Combinatorics

# ####################### MLE MVNormal ##############################
# d = 4
# bounds = easyBounds(d)
# nAgents = numAgents(d)
# agents = initializeAgents(nAgents, 1, d, bounds)
# C = .8
# F = .7
# # Define the mean vector and covariance matrix
# mean_vector = (1:d).^2
# covariance_matrix = Diagonal(ones(d))
# mv_normal = MvNormal(mean_vector, covariance_matrix)
# testFunc = x -> -1 * logpdf(mv_normal, vec(x))

# out = differentialEvolution(deepcopy(agents), C, F, testFunc, maxIters = 1000)
# print(out[1])
# print("\n")
# print(out[2])

# C = [.8 for x in 1:nAgents]
# F = [.7 for x in 1:nAgents]
# H = d * 2
# mF = [.7 for x in 1:H]
# mC = [.8 for x in 1:H]

# differentialEvolution(deepcopy(agents), C, F, testFunc, maxIters = 100, mF = mF, mC = mC, method = "current_to_pbest/1", variant = "SHADE", p = .2, H = H)
# differentialEvolution(deepcopy(agents), C, F, testFunc, maxIters = 100, mF = mF, mC = mC, method = "current_to_pbest/1", variant = "L-SHADE", p = .2, H = H)

# nArchive = Int(numAgents(d)/4)
# archive = initializeAgents(nArchive, 1, d, bounds)
# bounds = [[-1, 1] for i in 1:nAgents]

# differentialEvolution(deepcopy(agents), C, F, testFunc, maxIters = 100, mF = mF, mC = mC, method = "current_to_pbest/1", variant = "SHADE", p = .2, H = H, doArchive = true, archive = archive, nArchive = nArchive)
# differentialEvolution(deepcopy(agents), C, F, testFunc, maxIters = 100, mF = mF, mC = mC, method = "current_to_pbest/1", variant = "L-SHADE", p = .2, H = H, doArchive = true, archive = archive, nArchive = nArchive)

# differentialEvolution(deepcopy(agents), C, F, testFunc, maxIters = 100, mF = mF, mC = mC, method = "current_to_pbest/1", variant = "L-SHADE", p = .2, H = H, doArchive = true, archive = archive, nArchive = nArchive, bounds = bounds, boundType = "absorbingWall")

########################## D-optimal ###############################

K = 3
N = 16
d = K * N
bounds = easyBounds(K)
nAgents = 500
agents = initializeAgents(nAgents, N, K, bounds)
C = .8
F = .7
pArchive = 2.6
nArchive = Int(floor(nAgents * pArchive))
archive = initializeAgents(nArchive, N, K, bounds)
bounds = [[-1, 1] for i in 1:d]

# out = differentialEvolution(deepcopy(agents), C, F, testFunc, maxIters = 1000)
# print(out[1])
# print("\n")
# print(out[2])

C = [.8 for x in 1:nAgents]
F = [.7 for x in 1:nAgents]
H = 6
mF = [.5 for x in 1:H]
mC = [.5 for x in 1:H]

@time begin
out = differentialEvolution(deepcopy(agents), G_criterion, maxItersPSO = 1000, convergeIters = 1000, mF = mF, mC = mC, p = .11, H = H, bounds = bounds)
end

#45 seconds OG
#37 seconds with multithread
# print(out[1])
# print("\n")
# print(maxMinOptim(out[1]))
# plot(out[1][:, 1], out[1][:, 2], seriestype=:scatter)
# print("\n")
# print(out[2])
# print("\n")
# print(D_criterion(out[1]))
# print("\n")
# print(I_criterion(out[1]))
# print("\n")
# print(A_criterion(out[1]))
# print("\n")
print(G_criterion(out[1]))
print("\n")
print(out[4] / nAgents)