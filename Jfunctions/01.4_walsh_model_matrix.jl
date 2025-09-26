## =============================================================================
#
#   The purpose of this script is to take a design matrix as an input
#    (i.e. every candidate particle)
#
#

function genModelMat_fac(X; N, K ,order = 2)
    # NOTE: contrary to R's practice, I organanize the interaction terms before
    #       the squared terms === this is important for the I criterion
    #       wrt how the moments matrxi is organized.
    # X := design matrix
    # order : =  0 (first order no interactions)
    #         =  1 (first order interactions)
    #         =  2 (second order)
    #msize = size(X)
    #N     = msize[1]
    #k     = msize[2]

    # check order and k conditions


    # build model matrix
    if order == 0
        # grab columns and set model matrix for second order model
        if K == 1
            x1 = X[:, 1]
            Xm = [fill(1, N) x1]
        elseif K == 2
            x1 = X[:, 1]
            x2 = X[:, 2]
            Xm = [fill(1, N) x1 x2]
        elseif K == 3
            # min N = 10
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            Xm = [fill(1, N) x1 x2 x3]
        elseif K == 4
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            x4 = X[:, 4]
            Xm = [fill(1, N) x1 x2 x3 x4]
        elseif K == 5
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            x4 = X[:, 4]
            x5 = X[:, 5]
            Xm = [fill(1, N) x1 x2 x3 x4 x5]
        end
    elseif order == 1
        # grab columns and set model matrix for second order model
        if K == 1
            x1 = X[:, 1]
            Xm = [fill(1, N) x1]
        elseif K == 2
            x1 = X[:, 1]
            x2 = X[:, 2]
            Xm = [fill(1, N) x1 x2 (x1 .* x2)]
        elseif K == 3
            # min N = 10
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            Xm = [fill(1, N) x1 x2 x3 (x1 .* x2) (x1 .* x3) (x2 .* x3) ]
        elseif K == 4
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            x4 = X[:, 4]
            Xm = [fill(1, N) x1 x2 x3 x4 (x1 .* x2) (x1 .* x3) (x1 .* x4) (x2 .* x3) (x2 .* x4) (x3 .* x4)]
        elseif K == 5
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            x4 = X[:, 4]
            x5 = X[:, 5]
            Xm = [fill(1, N) x1 x2 x3 x4 x5 (x1 .* x2) (x1 .* x3) (x1 .* x4) (x1 .* x5) (x2 .* x3) (x2 .* x4) (x2 .* x5) (x3 .* x4) (x3 .* x5) (x4 .* x5)]
        end
    elseif order == 2
        # grab columns and set model matrix for second order model
        if K == 1
            x1 = X[:, 1]
            Xm = [fill(1, N) x1 x1.^2]
        elseif K == 2
            x1 = X[:, 1]
            x2 = X[:, 2]
            Xm = [fill(1, N) x1 x2 (x1 .* x2) x1.^2 x2.^2]
        elseif K == 3
            # min N = 10
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            Xm =  [fill(1, N) x1 x2 x3 (x1 .* x2) (x1 .* x3) (x2 .* x3) x1.^2 x2.^2 x3.^2]
        elseif K == 4
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            x4 = X[:, 4]
            Xm = [fill(1, N) x1 x2 x3 x4 (x1 .* x2) (x1 .* x3) (x1 .* x4) (x2 .* x3) (x2 .* x4) (x3 .* x4)  x1.^2 x2.^2 x3.^2 x4.^2]
        elseif K == 5
            x1 = X[:, 1]
            x2 = X[:, 2]
            x3 = X[:, 3]
            x4 = X[:, 4]
            x5 = X[:, 5]

            Xm = [fill(1, N) x1 x2 x3 x4 x5 (x1 .* x2) (x1 .* x3) (x1 .* x4) (x1 .* x5) (x2 .* x3) (x2 .* x4) (x2 .* x5) (x3 .* x4) (x3 .* x5) (x4 .* x5) x1.^2 x2.^2 x3.^2 x4.^2 x5.^2]
        end
    elseif order == 3
        # grab columns and set model matrix for second order model
        if K == 1
            x1 = X[:, 1]
            Xm = [fill(1, N) x1 x1.^2 x1.^3]
        elseif K == 2
            x1 = X[:, 1]
            x2 = X[:, 2]
            Xm = [fill(1, N) x1 x2 (x1 .* x2) x1.^2 x2.^2 (x1 .* x2.^2) x1.^3 x2.^3 (x1.^2 .* x2)]
        end
    end
    return Xm
end

# ## profile
#=
 using Distributions, Random, LinearAlgebra
Nt = 10
Kt = 2
Xt = genRandDesign_fac(; N = Nt, k = Kt, l_vec = fill(-1, Kt), u_vec = fill(1, Kt))
@time genModelMat_fac(Xt; N = Nt, K = Kt, order = 1)
=#
# # second run: 0.000038 seconds (29 allocations: 2.703 KiB)
