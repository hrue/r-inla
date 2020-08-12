checkit = function(n, eps = sqrt(.Machine$double.eps))
{
    Q.rw1.c = toeplitz(c(2, -1, rep(0, n-3), -1))
    Q.rw1 = Q.rw1.c
    Q.rw1[1, n] = Q.rw1[n, 1] = 0
    Q.rw1[1, 1] = Q.rw1[n, n] = 1

    Q.rw2.c = Q.rw1.c %*% t(Q.rw1.c)
    Q.rw2 = Q.rw1 %*% t(Q.rw1)
    Q.rw2[1, ] = c(1, -2, 1, rep(0, n-3))
    Q.rw2[2, ] = c(-2, 5, -4, 1, rep(0, n-4))
    Q.rw2[n, ] = rev(Q.rw2[1, ])
    Q.rw2[n-1, ] = rev(Q.rw2[2, ])

    A.rw1 = matrix(1, 1, n)
    A.rw1.c = matrix(1, 1, n)
    A.rw2 = rbind(rep(1, n), 1:n)
    A.rw2.c = matrix(1, 1, n)

    get.scale = function(Q, A, eps = sqrt(.Machine$double.eps))
    {
        rankdef = dim(A)[1]
        scale1 = exp(mean(log(diag(my.pinv(Q, rankdef=rankdef)))))

        if (is.null(A)) {
            Sig = solve(Q)
        } else {
            diag(Q) = diag(Q) + eps
            Qinv = solve(Q)
            Sig = Qinv - Qinv %*% t(A) %*% solve(A %*% Qinv %*% t(A)) %*% A %*% Qinv
        }
        scale2 = exp(mean(log(diag(Sig))))

        ## correction at the edges
        d = log(diag(Sig))
        n = dim(A)[2]
        scale3 = exp((0.5*d[1] + sum(d[2:(n-1)]) + 0.5*d[n])/(n-1))

        return(list(scale1=scale1, scale2=scale2, scale3=scale3))
    }

    return (list(
            rw1 = get.scale(Q.rw1, A.rw1, eps=eps), 
            rw1.c = get.scale(Q.rw1.c, A.rw1.c, eps=eps), 
            rw2 = get.scale(Q.rw2, A.rw2, eps=eps), 
            rw2.c = get.scale(Q.rw2.c, A.rw2.c, eps=eps)))
}


eps = .Machine$double.eps^0.75
scale.rw2.1 = c()
scale.rw2.2 = c()
for(n in seq(10, 400, by=20)) {
    print(n)
    a = checkit(n, eps=eps)
    scale.rw2.1[n] = a$rw2$scale1
    scale.rw2.2[n] = a$rw2$scale2
}
