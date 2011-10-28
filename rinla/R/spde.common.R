

inla.extract.el = function(M, ...)
{
    if (is.null(M))
        return(NULL)
    UseMethod("inla.extract.el", M)
}

inla.regex.match =  function(x, match) {
    return(strsplit(x, match)[[1]][1]=="")
}

inla.extract.el.matrix = function(M, match, by.row=TRUE)
{
    if (by.row) {
        return(M[sapply(rownames(M), inla.regex.match, match=match),,drop=FALSE])
    } else {
        return(M[,sapply(colnames(M), inla.regex.match, match=match),drop=FALSE])
    }
}

inla.extract.el.data.frame = function(M, match, by.row=TRUE)
{
    if (by.row) {
        return(M[sapply(rownames(M), inla.regex.match, match=match),,drop=FALSE])
    } else {
        return(M[,sapply(colnames(M), inla.regex.match, match=match),drop=FALSE])
    }
}

inla.extract.el.list = function(M, match)
{
    return(M[sapply(names(M), inla.regex.match, match=match)])
}



inla.spde.homogenise_B_matrix = function(B, n.spde, n.theta)
{
    if (!is.numeric(B))
        stop("B matrix must be numeric.")
    if (is.matrix(B)) {
        if ((nrow(B) != 1) && (nrow(B) != n.spde)) {
            stop(inla.paste(list("B matrix has",
                                 as.character(nrow(B)),
                                 "rows but should have 1 or",
                                 as.character(n.spde),
                                 sep=" ")))
        }
        if ((ncol(B) != 1) && (ncol(B) != 1+n.theta)) {
            stop(inla.paste(list("B matrix has",
                                 as.character(ncol(B)),
                                 "columns but should have 1 or",
                                 as.character(1+n.theta),
                                 sep=" ")))
        }
        if (ncol(B) == 1) {
            return(cbind(as.vector(B), matrix(0.0, n.spde, n.theta)))
        } else if (ncol(B) == 1+n.theta) {
            if (nrow(B) == 1) {
                return(matrix(as.vector(B), n.spde, 1+n.theta, byrow=TRUE))
            } else if (nrow(B) == n.spde) {
                return(B)
            }
        }
    } else { ## !is.matrix(B)
        if ((length(B) == 1) || (length(B) == n.spde)) {
            return(cbind(B, matrix(0.0, n.spde, n.theta)))
        } else if (length(B) == 1+n.theta) {
            return(matrix(B, n.spde, 1+n.theta, byrow=TRUE))
        } else {
            stop(inla.paste(list("Length of B vector is",
                                 as.character(length(B)),
                                 "but should be 1,",
                                 as.character(1+n.theta), "or",
                                 as.character(n.spde)),
                            sep=" "))
        }
    }
    stop(inla.paste(list("Unrecognised structure for B matrix"),
                    sep=" "))
}



inla.matern.cov = function(nu,kappa,x,d=1,corr=FALSE,theta, epsilon=1e-8)
{
    if (missing(theta)) { ## Ordinary Matern
        y = kappa*abs(x)
        if (corr) {
            ok = (y>=epsilon)
            if (nu<=0) {
                covariance = y*0
                covariance[!ok] = 1-y/epsilon
            } else {
                covariance = y*0
                covariance[ok] =
                    2^(1-nu)/gamma(nu) * (y[ok])^nu*besselK(y[ok], nu)
                if (any(!ok)) {
                    b = epsilon^(nu+1)*besselK(epsilon, nu-1)/
                        (gamma(nu)*2^(nu-1)-epsilon^nu*besselK(epsilon, nu))
                    yy = (y[!ok]/epsilon)^b
                    covariance[!ok] =
                        (1*(1-yy) +
                         (2^(1-nu)/gamma(nu)*
                          epsilon^nu*besselK(epsilon, nu))*yy)
                }
            }
            return(covariance)
        } else {
            ok = (y>=epsilon)
            covariance = y*0
            covariance[ok] =
                2^(1-nu)/gamma(nu+d/2)/(4*pi)^(d/2)/kappa^(2*nu)*
                    (y[ok])^nu*besselK(y[ok], nu)
            if (any(!ok)) {
                if (nu>0) { ## Regular Matern case
                    b = epsilon^(nu+1)*besselK(epsilon, nu-1)/
                        (gamma(nu)*2^(nu-1)-epsilon^nu*besselK(epsilon, nu))
                    yy = (y[!ok]/epsilon)^b
                    covariance[!ok] =
                        ((2^(1-nu)/gamma(nu+d/2)/(4*pi)^(d/2)/kappa^(2*nu))*
                         (gamma(nu)*2^(nu-1)*(1-yy) +
                          epsilon^nu*besselK(epsilon, nu)*yy))
                } else if (nu==0) { ## Limiting Matern case
                    g = 0.577215664901484 ## Euler's constant
                    covariance[!ok] =
                        2/gamma(d/2)/(4*pi)^(d/2)*
                            (-log(y[!ok]/2)-g)
                } else { ## (nu<0)
                    ## TODO: check this...
                    covariance[!ok] =
                        ((2^(1-nu)/gamma(nu+d/2)/(4*pi)^(d/2)/kappa^(2*nu)*
                          gamma(nu)*2^(nu-1))*(1-(y[!ok]/epsilon)) +
                         (2^(1-nu)/gamma(nu+d/2)/(4*pi)^(d/2)/kappa^(2*nu)*
                          epsilon^nu*besselK(epsilon, nu))*(y[!ok]/epsilon))
                }
            }
            return(covariance)
        }
    } else { ## Oscillating covariances
        y = abs(x)
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

        return(covariance*noise.variance)
    }
}


inla.matern.cov.s2 = function(nu,kappa,x,norm.corr=FALSE,theta=0)
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



inla.spde.models = function()
{
    types = c("spde1", "spde2")
    models = list()
    for (t in types) {
        models[[t]] =
            do.call(what=paste("inla.", t, ".models", sep=""),
                    args=list())
    }
    return(models)
}


inla.spde.sample =
    function(precision, seed=NULL)
{
    return(inla.finn(precision,
                     seed=(inla.ifelse(is.null(seed),
                                       0L,
                                       seed)))$sample)
}



inla.spde.precision = function(...)
{
    UseMethod("inla.spde.precision")
}

inla.spde.result = function(...)
{
    UseMethod("inla.spde.result")
}

