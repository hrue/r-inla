L = 1
n = 101
h = L/(n-1)
nE = 5
ntot = n+(nE-1)*2

if (nE>1) {
    c0E = nE*h*log((nE-((1:(nE-1))-0.5))/(nE-((1:(nE-1))+0.5)))
    c0 =
        c(c0E[(nE-1):1],
          h/2-nE*h*log((nE-0.5)/nE),
          rep(h, n-2),
          h/2-nE*h*log((nE-0.5)/nE),
          c0E[1:(nE-1)])
    gE = (1-2*(nE-(0:(nE-1))))/(2*nE*h)
    g1tmp =
        c(gE[nE:1],
          rep(-1/h, n-1),
          gE)
    g10 = -g1tmp[2:(ntot+1)]-g1tmp[1:ntot]
    g1p = g1tmp[2:ntot]

} else if (nE==1) {
    gE = -1/(2*nE*h)
    c0 =
        c(h/2-nE*h*log((nE-0.5)/nE),
          rep(h, n-2),
          h/2-nE*h*log((nE-0.5)/nE))
} else { ## nE==0
    c0 =
        c(h/2,
          rep(h, n-2),
          h/2)
}

C0 = Diagonal(ntot, c0)
G1 =
    sparseMatrix(i=c(1:ntot, 1:(ntot-1), 2:ntot),
                 j=c(1:ntot, 2:ntot, 1:(ntot-1)),
                 x=c(g10, g1p, g1p),
                 dims=c(ntot, ntot))


vari = diag(solve(C0*8/(L*0.2)^2+G1))
plot(((-nE+1):(n-1+(nE-1)))*h, vari)
lines(c(0,0), c(0,max(vari)*1.1))
lines(c(0,0)+L, c(0,max(vari)*1.1))
