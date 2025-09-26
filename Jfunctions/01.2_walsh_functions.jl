function D_criterion(X; N = N, K = K, order=2)
    Xm          = genModelMat_fac(X; N = N, K = K, order = order)
    msize       = size(Xm)
    #N           = msize[1]
    p           = msize[2]
    XpX         = transpose(Xm)*Xm
    determinant = det(XpX)

    if rank(XpX) < p
        result = -Inf
    else
        result = 100 * (determinant)^(1/p) / N
    end

    result = -result

    return result
end


## A criterion -----------------------------------------------------------------
function A_criterion(X; N = N, K = K,order=2)
    Xm          = genModelMat_fac(X; N = N, K = K, order = order)
    msize       = size(Xm)
    #N           = msize[1]
    p           = msize[2]
    XpX         = transpose(Xm)*Xm

    if rank(XpX) < p
        result = -Inf
    else
        XpX_inv = inv(XpX)
        result = 100*p/tr(N.*XpX_inv)
    end

    result = -result

    return result
end


## I criterion -----------------------------------------------------------------
function I_criterion(X; N = N, K = K, order=2)
    ## first moment matrices, then function eval
    # momen matrices are p x p
    # one factor moments matrix for second order model 3 x 3


    ## function eval
    Xm          = genModelMat_fac(X; N = N, K = K, order = order)
    msize       = size(Xm)
    #N           = msize[1]
    p           = msize[2]
    XpX         = transpose(Xm)*Xm
    # volume of design space 2^k
    #k = size(X)[2]
    V = 2^K
    # we need inv(XpX), first block matrices with small determinants
    determinant = det(XpX)

    if K == 1
        W = W1
    elseif K == 2
        W = W2
    elseif K == 3
        W =  W3
    elseif K == 4
        W = W4
    end

    if rank(XpX) < p
         result = -Inf
    else
         result = V / (N * tr(XpX \ W))
    end

    result = -result

    return result
end


## G-criterion -----------------------------------------------------------------
function G_criterion(X; N = N, K = K, order=2)
   ## function eval
    Xm          = genModelMat_fac(X; N = N, K = K, order = order)
    msize       = size(Xm)
    p           = msize[2]
    XpX         = transpose(Xm)*Xm

    if K == 1
        Xpred = Xs_onefac
    elseif K == 2
        Xpred = Xs_twofac
    elseif K == 3
        Xpred = Xs_threefac
    elseif K == 4
        Xpred = Xs_fourfac
    elseif K == 5
        Xpred = Xs_fivefac
        #println("size of Xpred grid is")
        #println(size(Xpred))
    end

    if rank(XpX) < p
         result = -Inf
    else
        #H = Xpred * (XpX \ transpose(Xpred))
        #H = Xpred * (\(XpX, transpose(Xpred)))
        #D = diag(N.*H)
        #Mx = maximum(D)
        C  = cholesky(XpX, check = false)
        Z  = \(transpose(C.U), transpose(Xpred))
        T  = @. N*Z^2
        D  = sum(T, dims = 1)
        Mx = maximum(D)
        result = 100*p/Mx
    end

    result = -result

    return result
end