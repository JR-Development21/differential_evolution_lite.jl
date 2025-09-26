## ================================================================================================
## function to make G grid

# Number of points (set at 11 in initial runs, JoBoGA uses 5)
Ng = 5
#Ng = 31
# JoBo suggested to sample the boundaries more densly, I used 21 initially
Nb = 5

## one factor -------------------------
x1 = collect(range(-1, stop = 1, length = Ng))
Xg = Matrix{Float64}(undef, length(x1), 3)
for i in 1:size(Xg)[1]
    Xg[i, :] = [1 x1[i] x1[i]^2]
end
Xs_onefac = Xg

## two factor ------------------------
x1  = collect(range(-1, stop = 1, length = Ng))
x2  = collect(range(-1, stop = 1, length = Ng))
Xg  = vec(collect(Iterators.product(x1, x2)))

Xg2 = Matrix{Float64}(undef, length(Xg), 2)
for i in 1:length(Xg)
    Xg2[i,:] = [Xg[i][1], Xg[i][2]]
end
#=
## make a finer grid on the boundary
Xg = vec(collect(Iterators.product([-1,1], collect(range(-1, stop = 1, length = Nb)))))
Xg3 = Matrix{Float64}(undef, length(Xg), 2)
for i in 1:length(Xg)
    Xg3[i,:] = [Xg[i][1], Xg[i][2]]
end
Xg = vec(collect(Iterators.product(collect(range(-1, stop = 1, length = Nb)),[-1,1])))
Xg4 = Matrix{Float64}(undef, length(Xg), 2)
for i in 1:length(Xg)
    Xg4[i,:] = [Xg[i][1], Xg[i][2]]
end
Xg = [Xg2; Xg3; Xg4]
Xg = [ones(size(Xg)[1]) Xg]
# remove duplicated rows
Xg = unique(Xg, dims = 1)
=#
Xg = [ones(size(Xg)[1]) Xg2]
x1 = Xg[:, 2]
x2 = Xg[:, 3]
Xg = [Xg x1.*x2 x1.^2 x2.^2]

Xs_twofac = Xg

## three fac grid
## k = 3 ---
x1  = collect(range(-1, stop = 1, length = Ng))
x2  = collect(range(-1, stop = 1, length = Ng))
x3  = collect(range(-1, stop = 1, length = Ng))

Xg  = vec(collect(Iterators.product(x1, x2, x3)))

Xg2 = Matrix{Float64}(undef, length(Xg), 3)
for i in 1:length(Xg)
    Xg2[i,:] = [Xg[i][1], Xg[i][2], Xg[i][3]]
end
#=
Xg = vec(collect(Iterators.product([-1,1],
    collect(range(-1, stop = 1, length = Nb)),
    collect(range(-1, stop = 1, length = Nb)))))

Xg3 = Matrix{Float64}(undef, length(Xg), 3)
for i in 1:length(Xg)
    Xg3[i,:] = [Xg[i][1], Xg[i][2], Xg[i][3]]
end
Xg =  vec(collect(Iterators.product(collect(range(-1, stop = 1, length = Nb)),
    collect(range(-1, stop = 1, length = Nb)),
    [-1,1])))
Xg4 = Matrix{Float64}(undef, length(Xg), 3)
for i in 1:length(Xg)
    Xg4[i,:] = [Xg[i][1], Xg[i][2], Xg[i][3]]
end
Xg =  vec(collect(Iterators.product(collect(range(-1, stop = 1, length = Nb)),
    [-1,1],
    collect(range(-1, stop = 1, length = Nb)))))
Xg5 = Matrix{Float64}(undef, length(Xg), 3)
for i in 1:length(Xg)
    Xg5[i,:] = [Xg[i][1], Xg[i][2], Xg[i][3]]
end

Xg = [Xg2; Xg3; Xg4; Xg5]
Xg = [ones(size(Xg)[1]) Xg]
# remove duplicated rows
Xg = unique(Xg, dims = 1)
=#
Xg = [ones(size(Xg)[1]) Xg2]
x1 = Xg[:, 2]
x2 = Xg[:, 3]
x3 = Xg[:, 4]
Xg = [Xg x1.*x2 x1.*x3 x2.*x3 x1.^2 x2.^2 x3.^2]
Xs_threefac = Xg


## four fac grid
## k = 4 ---
x1  = collect(range(-1, stop = 1, length = Ng))
x2  = collect(range(-1, stop = 1, length = Ng))
x3  = collect(range(-1, stop = 1, length = Ng))
x4  = collect(range(-1, stop = 1, length = Ng))

Xg  = vec(collect(Iterators.product(x1, x2, x3, x4)))
Xg2 = Matrix{Float64}(undef, length(Xg), 4)
for i in 1:length(Xg)
    Xg2[i,:] = [Xg[i][1], Xg[i][2], Xg[i][3], Xg[i][4]]
end

Xg = [ones(size(Xg)[1]) Xg2]

x1 = Xg[:, 2]
x2 = Xg[:, 3]
x3 = Xg[:, 4]
x4 = Xg[:, 5]
Xg = [Xg (x1 .* x2) (x1 .* x3) (x1 .* x4) (x2 .* x3) (x2 .* x4) (x3 .* x4) x1.^2 x2.^2 x3.^2 x4.^2]

Xs_fourfac = Xg

## five fac grid
## k = 5 ---
x1  = collect(range(-1, stop = 1, length = Ng))
x2  = collect(range(-1, stop = 1, length = Ng))
x3  = collect(range(-1, stop = 1, length = Ng))
x4  = collect(range(-1, stop = 1, length = Ng))
x5  = collect(range(-1, stop = 1, length = Ng))

Xg  = vec(collect(Iterators.product(x1, x2, x3, x4, x5)))
Xg2 = Matrix{Float64}(undef, length(Xg), 5)
for i in 1:length(Xg)
    Xg2[i,:] = [Xg[i][1], Xg[i][2], Xg[i][3], Xg[i][4], Xg[i][5]]
end

Xg = [ones(size(Xg)[1]) Xg2]

x1 = Xg[:, 2]
x2 = Xg[:, 3]
x3 = Xg[:, 4]
x4 = Xg[:, 5]
x5 = Xg[:, 6]

Xg = [Xg (x1 .* x2) (x1 .* x3) (x1 .* x4) (x1 .* x5) (x2 .* x3) (x2 .* x4) (x2 .* x5) (x3 .* x4) (x3 .* x5) (x4 .* x5) x1.^2 x2.^2 x3.^2 x4.^2 x5.^2]

Xs_fivefac = Xg


## I crit constants

W1 = [  2       0     2/3   ;
        0       2/3   0     ;
        2/3     0     2/5
    ]

# two factor moment matrix for second order model 6 x 6
W2 = [  4    0    0    0  4/3  4/3 ;
        0  4/3    0    0    0    0 ;
        0    0  4/3    0    0    0 ;
        0    0    0  4/9    0    0 ;
        4/3  0    0    0  4/5  4/9 ;
        4/3  0    0    0  4/9  4/5
    ]

 # three factor moment matrix for second order model 10 x 10
W3 = [  8   0   0   0   0   0   0   8/3 8/3 8/3 ;
        0   8/3 0   0   0   0   0   0   0   0   ;
        0   0   8/3 0   0   0   0   0   0   0   ;
        0   0   0   8/3 0   0   0   0   0   0   ;
        0   0   0   0   8/9 0   0   0   0   0   ;
        0   0   0   0   0   8/9 0   0   0   0   ;
        0   0   0   0   0   0   8/9 0   0   0   ;
        8/3 0   0   0   0   0   0   8/5 8/9 8/9 ;
        8/3 0   0   0   0   0   0   8/9 8/5 8/9 ;
        8/3 0   0   0   0   0   0   8/9 8/9 8/5
    ]

# four factor moment matrix for second order model 15 x 15
W4 = [  16    0    0    0    0    0    0    0    0    0    0 16/3 16/3 16/3 16/3 ;
        0 16/3    0    0    0    0    0    0    0    0    0    0    0    0    0 ;
        0    0 16/3    0    0    0    0    0    0    0    0    0    0    0    0 ;
        0    0    0 16/3    0    0    0    0    0    0    0    0    0    0    0 ;
        0    0    0    0 16/3    0    0    0    0    0    0    0    0    0    0 ;
        0    0    0    0    0 16/9    0    0    0    0    0    0    0    0    0 ;
        0    0    0    0    0    0 16/9    0    0    0    0    0    0    0    0 ;
        0    0    0    0    0    0    0 16/9    0    0    0    0    0    0    0 ;
        0    0    0    0    0    0    0    0 16/9    0    0    0    0    0    0 ;
        0    0    0    0    0    0    0    0    0 16/9    0    0    0    0    0 ;
        0    0    0    0    0    0    0    0    0    0 16/9    0    0    0    0 ;
     16/3    0    0    0    0    0    0    0    0    0    0 16/5 16/9 16/9 16/9 ;
     16/3    0    0    0    0    0    0    0    0    0    0 16/9 16/5 16/9 16/9 ;
     16/3    0    0    0    0    0    0    0    0    0    0 16/9 16/9 16/5 16/9 ;
     16/3    0    0    0    0    0    0    0    0    0    0 16/9 16/9 16/9 16/5  ]
