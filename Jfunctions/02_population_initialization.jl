############### Functions to initialize the population #####################

function numAgents(d; method = "linear")
    if method == "linear"
        return(10 * d)
    end
end

function initializeAgents(popSize, N, K, bounds; method = "uniform")
    if size(bounds)[3] != K
        throw("Mismatch between number of bounds given and the number of columns")   
    end
    agents = Array{Float64}(undef, (N, K, popSize))
    if method =="uniform"
        for i in 1:popSize
            for j in 1:K
                agents[:,j,i] = rand(Uniform(bounds[1,1,j], bounds[1,2,j]), N)
            end
        end
    end
    return agents
end

function easyBounds(nK)
    bounds = Array{Float64}(undef,(1,2,nK))
    bounds[:,1,:] .= -1
    bounds[:,2,:] .= 1
    return bounds
end