require(lattice)
require(Matrix)
require(orthopolynom)

plot.fake.points = function(s, radius=1, p=8, zoffset=0, ...)
{
    N = dim(s)[1]
    if (identical(dim(s)[2],2L)) {
      s = cbind(s,0)
    }

    offset = cbind(cos(2*pi*(0:(p-1))/p),sin(2*pi*(0:(p-1))/p),0)
    offset = t(matrix(rep(t(offset),N),3,p*N))*radius
    s.rep = matrix(rep(s,each=p),p*N,3) + offset

    idx = t(cbind(0,1:(p-2),2:(p-1)))
    idx = rep(idx,N)+rep((0:(N-1))*p+1,each=(p-2)*3)

    triangles3d(s.rep[idx,1],s.rep[idx,2],s.rep[idx,3]+zoffset, ...)
}

plot.rgl.map = function(map, z, border=FALSE, zoffset=0, ...)
{
    if (border) {
        rgl.lines(c(0,1),c(0,0),c(0,0))
        rgl.lines(c(1,1),c(0,1),c(0,0))
        rgl.lines(c(1,0),c(1,1),c(0,0))
        rgl.lines(c(0,0),c(1,0),c(0,0))
    }
    N = dim(map)[1]
    i = 1
    x = c()
    y = c()
    while(i < N)
    {
        n = map[i,2]
        x = c(x, NA, map[i + 1:n,1])
        y = c(y, NA, map[i + 1:n,2])
        i = i + n + 1
    }
    if (missing(z)) {
        z = rep(0,length(x))
        z[is.na(x)] = NA
    }
    rgl.linestrips(x,y,z+zoffset, ...)
}

plot.map = function(map, z, border=FALSE, ...)
{
#    if (border) {
#        rgl.lines(c(0,1),c(0,0),c(0,0))
#        rgl.lines(c(1,1),c(0,1),c(0,0))
#        rgl.lines(c(1,0),c(1,1),c(0,0))
#        rgl.lines(c(0,0),c(1,0),c(0,0))
#    }
    N = dim(map)[1]
    i = 1
    x = c()
    y = c()
    while(i < N)
    {
        n = map[i,2]
        x = c(x, NA, map[i + 1:n,1])
        y = c(y, NA, map[i + 1:n,2])
        i = i + n + 1
    }
    lines(x,y, ...)
}


calc.map = function(map, z, border=FALSE, ...)
{
#    if (border) {
#        rgl.lines(c(0,1),c(0,0),c(0,0))
#        rgl.lines(c(1,1),c(0,1),c(0,0))
#        rgl.lines(c(1,0),c(1,1),c(0,0))
#        rgl.lines(c(0,0),c(1,0),c(0,0))
#    }
    N = dim(map)[1]
    i = 1
    x = c()
    y = c()
    while(i < N)
    {
        n = map[i,2]
        x = c(x, NA, map[i + 1:n,1])
        y = c(y, NA, map[i + 1:n,2])
        i = i + n + 1
    }
    return (list(x=x,y=y))
}


points2mesh.calc = function(prefix,s)
{
    n = dim(s)[1]
    data = (inla.create.spde(prefix=prefix,
                             points2mesh=s))
    tv = (inla.read.fmesher.file(paste(prefix,
                                       "tv",
                                       sep="")))+1L
    ti = (inla.read.fmesher.file(paste(prefix,
                                      "p2m.t",
                                      sep="")))+1L
    ok = (ti > 0L)
    ti[ti == 0L] = NA
    b = (inla.read.fmesher.file(paste(prefix,
                                      "p2m.b",
                                      sep="")))
    ii = cbind(1:n,1:n,1:n)
    A = (sparseMatrix(dims=c(dim(s)[1],data$n),
                  i = as.vector(ii[ok,]),
                  j = as.vector(tv[ti[ok],]),
                  x = as.vector(b[ok,])))

    return (list(s=s,t=ti,b=b,A=A,ok=ok))
}

levelplotmap = function(..., mm) {
    panel.levelplot(...)
    panel.lines(mm$x, mm$y, col="black")
}



spdiag = function(a)
{
    b = as.vector(a)
    n = length(b)
    return(sparseMatrix(dims=c(n,n),i=1:n,j=1:n,x=b))
}




sparse.pattern.coarsen = function(A,factor)
{
    AA = as(A,"dgTMatrix")
    AAA = list(i=AA@i+1, j=AA@j+1,values=AA@x)

    sparseMatrix(i=ceiling(AAA$i/factor),j=ceiling(AAA$j/factor),x=1)
}



extract.rows = function(M,rowname)
{
    fn =  function(x,name) {
        return(strsplit(x,name,fixed=TRUE)[[1]][1]=="")
    }

    return(M[sapply(rownames(M),fn,name=rowname),,drop=FALSE])
}





make.plotdata = function(mapgrid,thedata)
{
    plotdata = as.vector(mapgrid$A %*% thedata)
    plotdata[!mapgrid$ok] = NA
    plotdata = (matrix(plotdata,length(mapgrid$az),length(mapgrid$el)))
    return(plotdata)
}


my.levelplot = function(mapgrid,plotdata,map=NULL,yrescale=2,color.palette,...)
{
    y.axis = c(-90,-60,-45,-30,-15,0,15,30,45,60,90)
    y.scale = list(at=sin(y.axis/180*pi)*yrescale,labels=y.axis)
    x.axis = (-4:4)/4*180
    x.scale = list(at=x.axis/180*pi,labels=x.axis)

    if (!missing(map)) {
        map$y = map$y*yrescale
    }
    yy = mapgrid$y*yrescale

    if (TRUE) {
        print(levelplot(row.values=mapgrid$x,column.values=yy,x=plotdata,
                        mm=list(x=map$x,y=map$y),
                        panel=levelplotmap,
                        col.regions=color.palette,
                        xlim=range(mapgrid$x),ylim=range(yy),aspect="iso",
                        contour=FALSE,cuts=101,labels=FALSE,pretty=TRUE,
                        xlab="Longitude",ylab="Latitude",
                        scales=list(y=y.scale,x=x.scale),...))
    }
##    print(contour(row.values=mapgrid$x,column.values=yy,x=plotdata))
}



matern.cov = function(nu,kappa,x,d=1,norm.corr=FALSE,theta)
{
    epsilon = 1e-8
    y = abs(x)

    if (missing(theta)) { ## Ordinary Matern
        covariance = 2^(1-nu)/gamma(nu+d/2)/(4*pi)^(d/2)/kappa^(2*nu)*
            (kappa*y+epsilon)^nu*besselK(kappa*y+epsilon,nu)

        if (norm.corr) {
            noise.variance = gamma(nu+d/2)/gamma(nu)*(4*pi)^(d/2)*kappa^(2*nu)
        } else {
            noise.variance = 1
        }
    } else { ## Oscillating covariances
        if (d>2L) {
            warning('Dimension > 2 not implemented for oscillating models.')
        }
        freq.max = 1000/max(y)
        freq.n = 10000
        w = seq(0,freq.max,length.out=freq.n)
        dw = w[2]-w[1]
        spec = 1/(2*pi)^d/(kappa^4+2*kappa^2*cos(pi*theta)*w^2+w^4)^((nu+d/2)/2)
        if (d==1L) {
            covariance = y*0+spec[1]*dw
        } else {
            covariance = y*0
        }
        for (k in 2:freq.n) {
            if (d==1L) {
                covariance = covariance+2*cos(y*w[k])*spec[k]*dw
            } else {
                covariance = covariance + w[k]*besselJ(y*w[k],0)*spec[k]*dw
            }
        }

        if (norm.corr) {
            noise.variance = 1/covariance[1]
        } else {
            noise.variance = 1
        }
    }

    return(covariance*noise.variance)
}


matern.cov.s2 = function(nu,kappa,x,norm.corr=FALSE,theta=0)
{
    y = cos(abs(x))

    freq.max = 40L
    freq.n = freq.max+1L
    w = 0L:freq.max
    spec = 1/(kappa^4+2*kappa^2*cos(pi*theta)*w*(w+1)+w^2*(w+1)^2)^((nu+1)/2)
    leg = legendre.polynomials(freq.max)
    covariance = y*0
    for (k in 1:freq.n) {
        covariance = (covariance + (2*w[k]+1)/(4*pi)*spec[k]*
                      polynomial.values(leg[k],y)[[1]])
    }

    if (norm.corr) {
        noise.variance = 1/covariance[1]
    } else {
        noise.variance = 1
    }

    return(covariance*noise.variance)
}



rbind.dgTMatrix = function(A,B) {
  if (class(A) != "dgTMatrix") {
    A = as(A,"dgTMatrix")
  }
  if (class(B) != "dgTMatrix") {
    B = as(B,"dgTMatrix")
  }
  return (sparseMatrix(i=c(A@i+1L,B@i+1L+A@Dim[1]),
                       j=c(A@j+1L,B@j+1L),
                       x=c(A@x,B@x),
                       dims=c(A@Dim[1]+B@Dim[1], max(A@Dim[2], B@Dim[2]))) )
}

rbind.dgCMatrix = function(A,B) {
  if (class(A) != "dgTMatrix") {
    A = as(A,"dgTMatrix")
  }
  if (class(B) != "dgTMatrix") {
    B = as(B,"dgTMatrix")
  }
  return (rbind.dgTMatrix(A,B))
}

cbind.dgTMatrix = function(A,B) {
  if (class(A) != "dgTMatrix") {
    A = as(A,"dgTMatrix")
  }
  if (class(B) != "dgTMatrix") {
    B = as(B,"dgTMatrix")
  }
  return (sparseMatrix(i=c(A@i+1L,B@i+1L),
                       j=c(A@j+1L,B@j+1L+A@Dim[2]),
                       x=c(A@x,B@x),
                       dims=c(max(A@Dim[1],B@Dim[1]), A@Dim[2]+B@Dim[2])) )
}

cbind.dgCMatrix = function(A,B) {
  if (class(A) != "dgTMatrix") {
    A = as(A,"dgTMatrix")
  }
  if (class(B) != "dgTMatrix") {
    B = as(B,"dgTMatrix")
  }
  return (cbind.dgTMatrix(A,B))
}

