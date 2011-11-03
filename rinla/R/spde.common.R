inla.dBind = function(...)
{
    A = list(...)
    if (length(A)<1)
        return(NULL)
    if (length(A)==1)
        return(A[[1]])
    B = A[[1]]
    for (k in 2:length(A)) {
        B = (rBind(cBind(B, Matrix(0, nrow(B), ncol(A[[k]]))),
                   cBind(Matrix(0, nrow(A[[k]]), ncol(B)), A[[k]])))
    }
    return(B)
}

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


inla.spde.sample = function(...)
{
    UseMethod("inla.spde.sample")
}

inla.spde.sample.default =
    function(precision, seed=NULL)
{
    return(inla.finn(precision,
                     seed=(inla.ifelse(is.null(seed),
                                       0L,
                                       seed)))$sample)
}

inla.spde.sample.inla.spde =
    function(spde, seed=NULL, ...)
{
    precision = inla.spde.precision(spde, ...)
    return(inla.spde.sample(precision, seed=seed))
}



inla.spde.precision = function(...)
{
    UseMethod("inla.spde.precision")
}

inla.spde.result = function(...)
{
    inla.require.inherits(list(...)[[1]], "inla", "First parameter")
    inla.require.inherits(list(...)[[2]], "character", "Second parameter")
    UseMethod("inla.spde.result", list(...)[[3]])
}







inla.spde.make.index = function(name, n.mesh, n.group=1, n.repl=1, n.field=n.mesh)
{
    if (!missing(n.field)) {
        warning("'n.field' is deprecated, please use 'n.mesh' instead.")
        if (missing(n.mesh))
            n.mesh = n.field
    }
    name.group = paste(name, ".group", sep="")
    name.repl = paste(name, ".repl", sep="")
    out = list()
    out[[name]]       = rep(rep(1:n.mesh, times=n.group), times=n.repl)
    out[[name.group]] = rep(rep(1:n.group, each=n.field), times=n.repl)
    out[[name.repl]]  = rep(1:n.repl, each=n.field*n.group)
    return(out)
}

inla.spde.make.A =
    function(mesh = NULL,
             loc = NULL,
             index = NULL,
             group = NULL,
             repl = 1L,
             n.mesh = NULL,
             n.group = max(group),
             n.repl = max(repl),
             group.mesh = NULL,
             group.method = c("nearest", "S0", "S1"))
{
    if (is.null(mesh)) {
        if (is.null(n.mesh))
            stop("At least one of 'mesh' and 'n.mesh' must be specified.")
    } else {
        inla.require.inherits(mesh, c("inla.mesh", "inla.mesh.1d"), "'mesh'")
        n.mesh = mesh$n
    }
    if (!is.null(group.mesh)) {
        inla.require.inherits(mesh, "inla.mesh.1d", "'mesh'")
    }
    group.method = match.arg(group.method)

    ## Handle loc and index input semantics:
    if (is.null(loc)) {
        A.loc = Diagonal(n.mesh, 1)
    } else {
        if (is.null(mesh))
            stop("'loc' specified but 'mesh' is NULL.")
        if (inherits(mesh, "inla.mesh.1d")) {
            A.loc = inla.mesh.1d.A(mesh, loc=loc, method="linear")
        } else {
            A.loc = inla.mesh.project(mesh, loc=loc)$A
        }
    }
    if (is.null(index)) {
        index = 1:nrow(A.loc)
    }
    ## Now 'index' points into the rows of 'A.loc'

    ## Handle group semantics:
    if (is.null(group.mesh)) {
        if (is.null(group))
            group = rep(1L, length(index))
        else if (length(group) == 1)
            group = rep(group, length(index))
        else if (length(group) != length(index))
            stop(paste("length(group) != length(index): ",
                       length(group), " != ", length(index),
                       sep=""))
    } else {
        n.group = group.mesh$n
        if (is.null(group))
            group = rep(mesh$loc[1], length(index))
        else if (length(group) == 1)
            group = rep(group, length(index))
        else if (length(group) != length(index))
            stop(paste("length(group) != length(index): ",
                       length(group), " != ", length(index),
                       sep=""))
        print(group)
        if (group.method=="nearest") {
            group.index =
                inla.mesh.1d.bary(group.mesh, loc=group, method="nearest")
            group = group.index$index[,1]
        } else {
            group.index =
                inla.mesh.1d.bary(group.mesh, loc=group, method="linear")
            if (group.method=="S0") {
                group = group.index$index[,1]
            }
        }
        print(group.index)
    }

    ## Handle repl semantics:
    if (is.null(repl))
        repl = rep(1, length(index))
    else if (length(repl) == 1)
        repl = rep(repl, length(index))
    else if (length(repl) != length(index))
        stop(paste("length(repl) != length(index): ",
                   length(repl), " != ", length(index),
                   sep=""))

    A.loc = inla.as.dgTMatrix(A.loc[index,,drop=FALSE])

    if (!is.null(group.mesh) && (group.method=="S1")) {
        return(sparseMatrix(i=(1L+c(A.loc@i, A.loc@i)),
                            j=(1L+c(A.loc@j+
                                    n.mesh*(group.index$index[,1]-1L)+
                                    n.mesh*n.group*(repl-1L),
                                    A.loc@j+
                                    n.mesh*(group.index$index[,2]-1L)+
                                    n.mesh*n.group*(repl-1L))),
                               x=c(A.loc@x*group.index$bary[,1],
                               A.loc@x*group.index$bary[,2]),
                            dims=c(length(index), n.mesh*n.group*n.repl)))
    } else {
        return(sparseMatrix(i=(1L+A.loc@i),
                            j=(1L+A.loc@j+
                               n.mesh*(group-1L)+
                               n.mesh*n.group*(repl-1L)),
                            x=A.loc@x,
                            dims=c(length(index), n.mesh*n.group*n.repl)))
    }
}


inla.stack = function(...)
{
    UseMethod("inla.stack")
}

inla.stack.default = function(data , A, effects, tag=NULL, strict=TRUE, ...)
{
    input.nrow = function(x) {
        return(inla.ifelse(is.matrix(x) || is(x, "Matrix"),
                           nrow(x),
                           length(x)))
    }
    input.ncol = function(x) {
        return(inla.ifelse(is.matrix(x) || is(x, "Matrix"),
                           ncol(x),
                           1))
    }
    input.list.ncol = function(l) {
        return(vapply(l,
                      function(x) inla.ifelse(is.list(x),
                                              input.ncol(x[[1]]),
                                              input.ncol(x)),
                      1))
    }
    input.list.nrow = function(l) {
        return(vapply(l,
                      function(x) inla.ifelse(is.list(x),
                                              input.nrow(x[[1]]),
                                              input.nrow(x)),
                      1))
    }
    effect.names = function(l) {
        islist = sapply(l, is.list)
        if (!is.null(names(l)))
            name = as.list(names(l)[!islist])
        else
            name = as.list(rep("", sum(!islist)))
        subnames =
            lapply(l,
                   function(x) inla.ifelse(is.list(x),
                                           names(x),
                                           NULL))
        return(setdiff(union(do.call(c, name),
                             do.call(c, subnames)),
                       ""))
    }
    expand.data.frames = function(l) {
        data.frames = vapply(l, is.data.frame, TRUE)
        return(c(l[!data.frames],
                 do.call(c, lapply(l[data.frames], as.list))))
    }
    effects.expand.data.frames = function(l) {
        result = l
        lists = vapply(l, is.list, TRUE)
        data.frames = vapply(l, is.data.frame, TRUE)
        for (k in which(!lists & !data.frames)) {
            result[[k]] = list(result[[k]])
            names(result[[k]])[1] = names(result)[k]
        }
        result[lists] = lapply(l[lists], expand.data.frames)
        result[data.frames] = lapply(l[data.frames], as.list)

        names(result) = rep("", length(result))
        return(result)
    }
    data.has.name = function(l) {
        if (is.null(names(l)))
            return(FALSE)
        return(all(vapply(names(l),
                          function(x) (x!=""),
                          TRUE)))
    }
    effect.has.name = function(l) {
        islist = sapply(l, is.list)
        if (!is.null(names(l)))
            name = as.list(names(l)[!islist])
        else
            name = as.list(rep("", sum(!islist)))
        subnames =
            lapply(l,
                   function(x) inla.ifelse(is.list(x),
                                           inla.ifelse(is.null(names(x)),
                                                       rep("", length(x)),
                                                       names(x)),
                                           NULL))
        return(all(vapply(c(do.call(c, name), do.call(c, subnames)),
                          function(x) (x!=""),
                          TRUE)))
    }
    effect.nrow = function(l) {
        return(vapply(l,
                      function(x) inla.ifelse(is.list(x),
                                              max(input.list.nrow(x)),
                                              input.nrow(x)),
                      1))
    }

    if (length(list(...))>0)
        warning(paste("Extra argument '", names(list(...)), "' ignored.",
                      collapse="\n", sep=""))

    ## Check if only a single block was specified.
    if (!is.list(A)) {
        A = list(A)
        effects =
            inla.ifelse(is.data.frame(effects),
                        list(expand.data.frames(list(effects))),
                        effects.expand.data.frames(list(effects)))
    } else {
        effects = effects.expand.data.frames(effects)
    }
    data = expand.data.frames(data)

    if (length(A) != length(effects))
        stop(paste("length(A)=", length(A),
                   " should be equal to length(effects)=", length(effects), sep=""))

    if (!data.has.name(data)) {
        stop(paste("All data must have names\n",
                   "Data names:   ",
                   paste(names(data), collapse=", ", sep=""),
                   "\n",
                   "Effect names: ",
                   paste(effect.names(effects), collapse=", ", sep=""),
                   sep=""))
    }
    if (!effect.has.name(effects)) {
        stop(paste("All effects must have names\n",
                   "Data names:   ",
                   paste(names(data), collapse=", ", sep=""),
                   "\n",
                   "Effect names: ",
                   paste(effect.names(effects), collapse=", ", sep=""),
                   sep=""))
    }
    if (length(unique(c(names(data), effect.names(effects)))) <
        length(data)+length(effect.names(effects))) {
        stop(paste("Names for data and effects must not coincide.\n",
                   "Data names:   ",
                   paste(names(data), collapse=", ", sep=""),
                   "\n",
                   "Effect names: ",
                   paste(effect.names(effects), collapse=", ", sep=""),
                   sep=""))
    }

    n = effect.nrow(effects)

    A.n = length(A)
    A.ncol = input.list.ncol(A)
    A.nrow = input.list.nrow(A)
    if (any((A.ncol==1) & (A.nrow==1) & (n!=1))) {
        idx = which((A.ncol==1) & (A.nrow==1) & (n!=1))
        for (i in idx) {
            A[[i]] = Diagonal(n[i], A[[i]])
        }
    }
    A.ncol = input.list.ncol(A)
    A.nrow = input.list.nrow(A)
    if (any((A.ncol==1) & (A.nrow==1) & (n==1))) {
        idx = which((A.ncol==1) & (A.nrow==1) & (n==1))
        for (i in idx) {
            A[[i]] = Matrix(A[[i]], max(A.nrow), 1)
        }
    }
    A.ncol = input.list.ncol(A)
    A.nrow = input.list.nrow(A)
    if (any(A.nrow < max(A.nrow))) {
        stop(paste("Mismatching nrow(A) for A[[",
                   paste(which(A.nrow != max(A.nrow)), collapse=","), "]]: ",
                   paste(A.nrow[A.nrow != max(A.nrow)], collapse=","),
                   " should be ", max(A.nrow), sep=""))
    }
    A.nrow = max(A.nrow)

    if (any(n!=A.ncol)) {
        stop(paste("Mismatching ncol(A) vs nrow(effects) for ",
                   paste(which(n!=A.ncol), collapse=","), ": ",
                   paste(A.ncol[n!=A.ncol], collapse=","),
                   " should be ",
                   paste(n[n!=A.ncol], collapse=","),
                   sep=""))
    }

    data.nrow = input.list.nrow(data)
    data.ncol = input.list.ncol(data)
    if (any(data.nrow > A.nrow)) {
        if (strict) {
            stop(paste("Strict mode:\n",
                       "Length of data '",
                       names(data[data.nrow > A.nrow]),
                       "' is ",
                       data.nrow[data.nrow > A.nrow],
                       " but should be ", A.nrow,
                       collapse="\n",
                       sep=""))
        } else {
            stop(paste("Length of data '",
                       names(data[data.nrow > A.nrow]),
                       "' is ",
                       data.nrow[data.nrow > A.nrow],
                       " but should be <=", A.nrow,
                       collapse="\n",
                       sep=""))
        }
    }

    data.output = data
    for (k in 1:length(data)) {
        if (data.nrow[k] < A.nrow) {
            if (strict) {
                stop(paste("Strict mode:\n",
                           "Length of data '",
                           names(data[data.nrow < A.nrow]),
                           "' is ",
                           data.nrow[data.nrow < A.nrow],
                           " but should be ", A.nrow,
                           collapse="\n",
                           sep=""))
            }
            if (is.matrix(data[[k]])) {
                data.output[[k]] =
                    rbind(data[[k]],
                          matrix(NA, A.nrow-data.nrow[k], data.ncol[k]))
            } else {
                data.output[[k]] =
                    c(data[[k]],
                      rep(NA, A.nrow-data.nrow[k]))
            }
        }
    }

    effects.output = list()
    n0 = 1
    ntot = sum(n)
    for (i in 1:length(effects)) {
        n1 = n0+n[i]-1
        nj = input.list.nrow(effects[[i]])
        if (strict && (any(nj != n[i]))) {
            stop("Strict mode:\n",
                 "Length of effect '", names(effects[[i]])[nj != n[i]],
                 "' is ",
                 nj[nj != n[i]],
                 " but should be ",
                 n[i], ".",
                 collapse = "\n",
                 sep="")
        }
        nmissing = n[i]-nj
        for (j in 1:length(effects[[i]])) {
            thename = names(effects[[i]])[j]
            if (is.matrix(effects[[i]][[j]])) {
                n.col = ncol(effects[[i]][[j]])
                if (is.null(effects.output[[thename]])) {
                    effects.output[[thename]] =
                    rbind(matrix(NA, n0-1, n.col),
                          effects[[i]][[j]],
                          matrix(NA, nmissing[j]+ntot-n1, n.col))
                } else {
                    if (!is.matrix(effects.output[[thename]]))
                        stop(paste("Effect '", thename,
                                   "' has inconsistent matrix property.",
                                   sep=""))
                    if (n.col != ncol(effects.output[[thename]]))
                        stop(paste("Effect '", thename,
                                   "' has inconsistent matrix columns: ",
                                   n.col,
                                   " previously specified as ",
                                   ncol(effects.output[[thename]]), ".",
                                   sep=""))
                    effects.output[[thename]][n0:n1,] =
                        rbind(effects[[i]][[j]],
                              matrix(NA, nmissing[j], n.col))
                }
            } else {
                if (is.null(effects.output[[thename]])) {
                    effects.output[[thename]] =
                        c(rep(NA, n0-1),
                          effects[[i]][[j]],
                          rep(NA, nmissing[j]+ntot-n1))
                } else {
                    effects.output[[thename]][n0:n1] =
                        c(effects[[i]][[j]],
                          rep(NA, nmissing[j]))
                }
            }
        }
        n0 = n1+1
    }

    A.output = do.call(cBind, A)

    index = list(1:A.nrow)
    if (!is.null(tag)) {
        names(index) = tag
    }

    stack =
        list(data=c(data.output, effects.output), A=A.output,
             data.names = names(data.output),
             effect.names = names(effects.output),
             n.data = A.nrow,
             index = index)
    class(stack) = "inla.data.stack"

    return(stack)
}

inla.stack.inla.data.stack = function(...)
{
    input.nrow = function(x) {
        return(inla.ifelse(is.matrix(x) || is(x, "Matrix"),
                           nrow(x),
                           length(x)))
    }
    input.ncol = function(x) {
        return(inla.ifelse(is.matrix(x) || is(x, "Matrix"),
                           ncol(x),
                           1))
    }
    input.list.ncol = function(l) {
        return(sapply(l,
                      function(x) inla.ifelse(is.list(x),
                                              input.ncol(x[[1]]),
                                              input.ncol(x))))
    }
    input.list.nrow = function(l) {
        return(sapply(l,
                      function(x) inla.ifelse(is.list(x),
                                              input.nrow(x[[1]]),
                                              input.nrow(x))))
    }
    input.parse = function(l) {
        return(sapply(l,
                      function(x) inla.ifelse(is.list(x),
                                              input.nrow(x[[1]]),
                                              input.nrow(x))))
    }

    S = list(...)
    if (length(S)<1)
        return(NULL)
    inla.require.inherits(S[[1]], "inla.data.stack", "Argument 1")
    if (length(S)==1)
        return(S[[1]])

    S1 = S[[1]]

    for (k in 2:length(S)) {
        inla.require.inherits(S[[k]], "inla.data.stack",
                              paste("Argument ", k, sep=""))


        ## rbind all the elements of S1 and S[[k]]
        all.data.names = union(S1$data.names, S[[k]]$data.names)
        all.effect.names = union(S1$effect.names, S[[k]]$effect.names)
        if (length(intersect(all.data.names, all.effect.names))>0) {
            stop("Name conflict.")
        }
        all.names = c(all.data.names, all.effect.names)

        present1 =
            is.element(all.names, c(S1$data.names, S1$effect.names))
        present2 =
            is.element(all.names, c(S[[k]]$data.names, S[[k]]$effect.names))
        ismat1 =
            lapply(S1$data[c(S1$data.names, S1$effect.names)],
                   is.matrix)
        ismat2 =
            lapply(S[[k]]$data[c(S[[k]]$data.names, S[[k]]$effect.names)],
                   is.matrix)

        data = list()
        for (j in 1:length(all.names)) {
            name = all.names[j]
            is.data = (is.element(name, S1$data.names) ||
                       is.element(name, S[[k]]$data.names))
            if (is.data) {
                size1 = S1$n.data
                size2 = S[[k]]$n.data
            } else {
                size1 = ncol(S1$A)
                size2 = ncol(S[[k]]$A)
            }
            if (present1[j] && present2[j]) {
                data[[name]] =
                    inla.ifelse(ismat1[[name]] || ismat2[[name]],
                                rbind(inla.ifelse(ismat1[[name]],
                                                  S1$data[[name]],
                                                  as.matrix(S1$data[[name]])),
                                      inla.ifelse(ismat2[[name]],
                                                  S[[k]]$data[[name]],
                                                  as.matrix(S[[k]]$data[[name]]))),
                                c(S1$data[[name]],
                                  S[[k]]$data[[name]]))
            } else if (present1[j]) {
                data[[name]] =
                    inla.ifelse(ismat1[[name]],
                                rbind(S1$data[[name]],
                                      matrix(NA, size2,
                                             ncol(S1$data[[name]]))),
                                c(S1$data[[name]],
                                  rep(NA, size2)))
            } else {
                data[[name]] =
                    inla.ifelse(ismat2[[name]],
                                rbind(matrix(NA, size1,
                                             ncol(S[[k]]$data[[name]])),
                                      S[[k]]$data[[name]]),
                                c(rep(NA, size1),
                                      S[[k]]$data[[name]]))
            }
##            if ((present1[j] && ismat1[[name]]) ||
##                (present2[j] && ismat2[[name]]))
##                print(dim(data[[name]]))
##            else
##                print(length(data[[name]]))
        }

        A = inla.dBind(S1$A, S[[k]]$A)
        index =
            c(S1$index,
              lapply(S[[k]]$index,
                     function(x) (x+S1$n.data)))
        n.data = S1$n.data+S[[k]]$n.data

        S1 =
            list(data=data, A=A,
                 data.names = all.data.names,
                 effect.names = all.effect.names,
                 n.data = S1$n.data + S[[k]]$n.data,
                 index = index)
        class(S1) = "inla.data.stack"
    }

    return(S1)
}



inla.stack.data = function(stack, ...)
{
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    return(c(stack$data, list(...)))
}

inla.stack.A = function(stack)
{
    inla.require.inherits(stack, "inla.data.stack", "'stack'")
    return(stack$A)
}
