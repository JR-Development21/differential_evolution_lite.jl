function D_optim(mat; order = 2)
    size_mat = size(mat)
    if length(size_mat) > 1
        K = size_mat[2]
    else
        K = 1
    end
    N = size_mat[1]
    X = createModelMatrix(mat, order, N, K)
    result = -det(X' * X)
    return result
end

function safe_func(func, A, agentLength, C, F)
    try
        return func(A)
    catch e
        println("Error computing rank: ", e)
        println("Matrix size: ", size(A))
        println("Contains NaN: ", any(isnan, A))
        println("agentLength: ", agentLength)
        println("F: ", F)
        println("C: ", C)
        println(A)
        println("Contains Inf: ", any(isinf, A))
        println("Condition number: ", cond(A))
        println("Singular values (if computable): ")
        try
            println(svdvals(A))
        catch svd_err
            println("Failed to compute singular values: ", svd_err)
        end
        return missing  # or return a default value
    end
end

function I_optim(mat; order = 2)
    size_mat = size(mat)
    if length(size_mat) > 1
        K = size_mat[2]
    else
        K = 1
    end
    N = size_mat[1]
    X = createModelMatrix(mat, order, N, K)
    V = 2^K
    W = genW(K)
    XX = X'X
    if det(XX) < sqrt(eps())
        result = Inf
    else
        result = try
            N * tr(W / XX) / V
        catch e
            Inf
        end
    end
    return(result)
end

function A_optim(mat; order = 2)
    size_mat = size(mat)
    if length(size_mat) > 1
        K = size_mat[2]
    else
        K = 1
    end
    N = size_mat[1]
    X = createModelMatrix(mat, order, N, K)
    XX = X'X
    if det(XX) < sqrt(eps())
        result = Inf
    else
        result = try
            tr(inv(XX))
        catch e
            Inf
        end
    end
    return(result)
end

function G_optim(mat; order = 2)
    size_mat = size(mat)
    if length(size_mat) > 1
        K = size_mat[2]
    else
        K = 1
    end
    N = size_mat[1]
    X = createModelMatrix(mat, order, N, K)
    XX = X'X
    if det(XX) < sqrt(eps())
        result = Inf
    else
        X_inv = try
            inv(XX)
        catch e
            "inverse failed"
        end
        if typeof(X_inv) == String
            result = Inf
        else
            base_vec = [[-1, -.5, 0, .5, 1] for i in 1:K]
            base_tuples = vec(collect(Base.product(base_vec...)))
            base_grid = covertToTuple(base_tuples)
            base_grid = createModelMatrix(base_grid, order, 5^K, K)
            max_track = -Inf
            for i in 1:(5^K)
                G_cand = N * base_grid[i,:]' * X_inv * base_grid[i,:]
                if G_cand > max_track
                    max_track = copy(G_cand)
                end
            end
            result = max_track
        end
    end
    return result
end

function maxMinOptim(mat)
    D = pairwise(Euclidean(), mat, dims=1)
    for i in axes(D, 1)
        D[i,i] = Inf
    end

    result = -minimum(D)
    return result
end

function genW(K)
    if K == 1
        M = 2 .* [1 0 1/3;
        0 1/3 0;
        1/3 0 1/5]
    else
        M = 2^K .* [1 zeros(K)' zeros(Int(K * (K - 1) / 2))' 1/3*ones(K)';
        zeros(K) Diagonal(1/3 * ones(K)) zeros(Int, K, Int(K * (K - 1) / 2)) zeros(Int, K, K);
        zeros(Int(K * (K - 1) / 2)) zeros(Int, Int(K * (K - 1) / 2), K) Diagonal(1/9 * ones(Int(K * (K-1) / 2))) zeros(Int, Int(K * (K - 1) / 2), K);
        1/3*ones(K) zeros(Int, K, K) zeros(Int, K, Int(K * (K - 1) / 2)) (Diagonal(1/5 * ones(K)) + 1/9 * (ones(K,K) - Diagonal(ones(K))))]
    end
    return M
end

function covertToTuple(x)
    out = Matrix{Float64}(undef, length(x), length(x[1]))
    for i in 1:length(x)
          for j in 1:length(x[1])
              out[i, j] = x[i][j]
          end
    end
    out
 end


function createModelMatrix(X, order::Int, N, K)

    if order < 1 || order > 3
        throw(ArgumentError("Only supports models of order 1 (linear), 2 (quadratic), and 3 (cubic)."))
    end

    # Start with an intercept column of 1's
    model_matrix = ones(N, 1)

    # Add linear terms (main effects)
    if order >= 1
        model_matrix = hcat(model_matrix, X)
    end

    # Add interaction terms and quadratic terms
    if order >= 2
        # Add interaction terms (pairwise interactions)
        for j in 1:K-1
            for k in j+1:K
                model_matrix = hcat(model_matrix, X[:, j] .* X[:, k])
            end
        end

        # Add quadratic terms (squared terms)
        for j in 1:K
            model_matrix = hcat(model_matrix, X[:, j] .^ 2)
        end
    end

    # Add cubic terms
    if order == 3
        for j in 1:K
            model_matrix = hcat(model_matrix, X[:, j] .^ 3)
        end

        # Add interaction of quadratic with linear (X_i^2 * X_j)
        for j in 1:K-1
            for k in j+1:K
                model_matrix = hcat(model_matrix, (X[:, j] .^ 2) .* X[:, k])
                model_matrix = hcat(model_matrix, X[:, j] .* (X[:, k] .^ 2))
            end
        end
    end

    return model_matrix
end

function DToBork(val, N, P)
    (-val)^(1/P) * 100 / N
end

function IToBork(val)
    1 / val
end

function AToBork(val, N, P)
    100 * P / val / N
end

function GToBork(val, P)
    100 * P / val
end

function design2Space(designVec, minVal, maxVal)
    (designVec .+ 1) ./ 2  .* (maxVal - minVal).+ minVal
end

# using LinearAlgebra

# test1 = [1;-1;0]

# (-D_optim(test1))^(1/3) * 100 / 3
# 1 / I_optim(test1)
# 100 * 3 / A_optim(test1) / 3
# 100 * 3 / G_optim(test1)

# test2 = [1;-1;0;0]

# print((-D_optim(test2))^(1/3) * 100 / 4)
# print(1 / I_optim(test2))
# print(100 * 3 / A_optim(test2) / 4)

# test3 = [1;-1;.4858682;-.4858682]

# print(100 * 3 / G_optim(test3))

# test4 = [1 1; -1 1; -1 -1; 1 -.394449; .394449 -1; -.131483 .131483]
# print((-D_optim(test4))^(1/6) * 100 / 6)

# test5 = [1 -1; -1 -1; 0 1; -1 .503816 ; 1 .503816 ; 0 -.220484]
# 100 * 6 / A_optim(test5) / 6

# test6 = [-1 -1; -.193256 -.193256; 1 .522; .522 1; .864001 -1; -1 .864001]
# 100 * 6 / G_optim(test6)

# test7 = [-1 1; .707479 1; -1 -.707479; 1 -.276367; .276367 -1; -.144868 .144868]
# 1 / I_optim(test7)

# test8 = [.2912 -1 -1; -1 .2912 -1; -1 -1 .2912; -.1925 1 -.1925; -.1925 -.1925 1; 1 -.1925 -.1925; -1 1 1; 1 1 -1; 1 -1 1; 1 1 1]
# print((-D_optim(test8))^(1/10) * 100 / 10)

# test9 = [-.1749 -1 1; 1 .1749 1; 1 -1 -.1749; -1 -.1072 .1072 ; .1072 1 .1072 ; .1072 -.1072 -1; -1 1 -1; -1 1 1; 1 1 -1; -1 -1 -1]
# 100 * 10 / A_optim(test9) / 10

# test10 = [1 0 0; 0 1 0; 0 0 1; 1 1 1; -1 1 1; 1 -1 1; 1 1 -1; -1 -1 1; -1 1 -1; 1 -1 -1; -1 -1 -1]
# print((-D_optim(test10))^(1/10) * 100 / 11)

# test11 = [.1223 -1 1; -1 .1223 1; -1 -1 -.1223; -.2166  -.2166  -1 ; -.2166  1  .2166  ; 1 -.2166  .2166 ; .0914 .0914 -.0914; 1 -1 -1; 1 1 -1; -1 1 -1; 1 1 1]
# 100 * 10 / A_optim(test11) / 11

# test12 = [-.1140 -1 1; -1 -.1140 1; -1 -1 .1140; -.7550 -.7550 -1; -.7550 1 .7550; 1 -.7550 .7550; .1858 .1858 -.1858; 1 -1 -1; 1 1 -1; -1 1 -1; 1 1 1]
# 100 * 10 / G_optim(test12)

# test13 = [-.0589 .0589 -.0589; -1 -.1905 .1905;.1905 1 .1905; .1905 -.1905 -1; -.2349 -1 1; 1 .2349 1; 1 -1 -.2349; -1 -1 -1; -1 1 -1; -1 1 1; 1 1 -1]
# 1 / I_optim(test13)
