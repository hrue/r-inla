## Export: inla.knmodels

##! \name{inla.knmodels}
##! \alias{inla.knmodels}
##! \alias{knmodels}
##! \title{Spacetime interaction models}
##! \description{
##!    It implements the models in Knorr-Held, L. (2000) 
##!    with three different constraint approaches: 
##!    sum-to-zero, contrast or diagonal add.
##! } 
##! \usage{
##!inla.knmodels(
##! formula,
##! data,
##! progress=FALSE, 
##! control.st=list(
##!   t=NULL, 
##!   s=NULL,
##!   st=NULL,
##!   graph=NULL,
##!   type=c(paste(1:4), paste0(2:4, 'c'), paste0(2:4, 'd')), 
##!   diagonal=1e-5, 
##!   ...) 
##! )
##!}
##! \arguments{
`inla.knmodels` =
    function(
       ##! \item{formula}{The formula specifying the other 
       ##!   model components, without the spacetime 
       ##!   interaction term. The spacetime interaction term 
       ##!   will be added accordly to the specification in 
       ##!   the \code{control.st} argument. See \code{inla}}
       formula, 
       ##! \item{progress}{If it is to be shown the model 
       ##!   fitting progress. Useful if more than one 
       ##!   interaction type is being fitted.}
       progress=FALSE,
       ##! \item{control.st}{Named list of arguments to control
       ##!   the spacetime interaction. It should contains: 
       control.st=list(
           ##!  \code{time} to be used as the index set for the
           ##!   main temporal effect which will be considered
           ##!   for the constraints when it is the case.
           time,
           ##!  \code{space} to be used as the index set for the
           ##!   main spatial effect which will be considered
           ##!   for the constraints when it is the case.
           space,
           ##!  \code{spacetime} to be the index set for the
           ##!   spacetime interaction effect.
           spacetime,
           ##! \code{graph} to be the graph for the spatial neighbor 
           ##!   structure to be used in a \code{\link{f}} term 
           ##!   for the main spatial random effect term or for 
           ##!   building the spacetime interaction model.
           graph,
           ##! \code{type} to specify the spacetime interaction type.  
           ##!   \code{1} to \code{4} corresponds to the four 
           ##!   interaction types in Knorr-Held, L. (2000) with 
           ##!   all the needed sum-to-zero constraints. 
           ##!   \code{2c}, \code{3c} and \code{4c} are 
           ##!   the contrast version considering the first time  
           ##!   or space constrained to be equal to zero. 
           ##!   \code{2d}, \code{3d} and \code{4d} are the 
           ##!   corresponding versions when considering the 
           ##!   diagonal add approach.  
           type=c(paste(1:4), paste0(2:4, 'c'), paste0(2:4, 'd')), 
           ##! \code{diagonal} to be the value to be added to the 
           ##!   diagonal when using the diagonal add approach.
           diagonal=1e-5, 
           ##! \code{timeref} to specify the time point to be the
           ##!   reference time in the contrast parametrization.
           timeref=1, 
           ##! \item{spaceref} to specify the area to be the
           ##!   reference for the contrast parametrization.
           spaceref=1, 
           ##!  \code{...} where additional arguments can be 
           ##!   passed to \code{\link{f}} function.
           ##!   Specification of the hyperparameter, 
           ##!   fixed or random, initial value, prior and its 
           ##!   parameters for the spacetime interaction. See 
           ##!   \code{?inla.models} and look for \code{generic0}. 
           ##!   By default we scale it and use the PC-prior to set 
           ##!   the prior using the \code{pc.prec} prior with 
           ##!   \code{param = c(0.5, 0.5)}. See documentation with 
           ##!   \code{?inla.doc("pc.prec")}.
           ...),
       ##! \item{...}{Arguments to be passed to the 
       ##!   \code{\link{inla}} function.}
       ...)
{
##! }
##!}
##! \value{
##!  \code{inla.knmodels} returns an object of class \code{"inla"}. 
##!    or a list of objects of this class if it is asked to compute 
##!    more than one interaction type at once. 
##! Note: when the model type is 2c, 3c, 4c, 2d, 3d or 4d, it also 
##!   includes linear combinations summary.
##! }
##! \author{Elias T. Krainski}
##! \seealso{
##!     \code{\link{inla.knmodels.sample}} to sample from
##! }
##! \examples{
##!### define space domain as a grid
##!grid <- SpatialGrid(GridTopology(c(0,0), c(1, 1), c(4, 5)))
##!(n <- nrow(xy <- coordinates(grid)))
##!
##!### build a spatial neighborhood list
##!jj <- lapply(1:n, function(i) 
##!    which(sqrt((xy[i,1]-xy[,1])^2 + (xy[i,2]-xy[,2])^2)==1))
##!
##!### build the spatial adjacency matrix
##!graph <- sparseMatrix(rep(1:n, sapply(jj, length)),
##!                      unlist(jj), x=1, dims=c(n, n))
##!
##!### some random data at 10 time points
##!dat <- inla.knmodels.sample(graph, m=10, tau.t=2, tau.s=2, tau.st=3)
##!str(dat)
##!sapply(dat$x, summary)
##!
##!nd <- length(dat$x$eta)
##!dat$e <- runif(nd, 0.9, 1.1)*rgamma(n, 40, 2)
##!dat$y <- rpois(nd, dat$e*exp(dat$x$eta-3))
##!summary(dat$y)
##!
##!### fit the type 4 considering three different approaches 
##!tgraph <- sparseMatrix(i=c(2:10, 1:9), j=c(1:9, 2:10), x=-1)
##!res <- inla.knmodels(y ~ f(time, model='bym2', graph=tgraph) +
##!     f(space, model='bym2', graph=graph),
##!     data=dat, family='poisson', E=dat$E, progress=TRUE, 
##!     control.st=list(time=time, space=space, 
##!        spacetime=spacetime, graph=graph, type=c(4, '4c', '4d')), 
##!     control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE))
##!sapply(res, function(x)
##!       c(dic=x$dic$dic, waic=x$waic$waic, cpo=-sum(log(x$cpo$cpo))))
##!}
    mcall <- match.call(expand.dots=TRUE)
    ft <- paste('~', mcall$control.st$time)
    if (ft=='~ ') {
        time <- tname <- NULL
    } else {
        tname <- substring(ft, 3)
        time <- model.frame(as.formula(ft), data=eval(mcall$data))[,1]
        timeref <- unique(eval(mcall$control.st$timeref))
        if (length(timeref)>1)
            stop("length(timeref)>1")
        if (length(timeref)==0) timeref <- 1
    }
    fs <- paste('~', mcall$control.st$space)
    if (fs=='~ ') {
        space <- sname <- NULL
    } else {
        sname <- substring(fs, 3)
        space <- model.frame(as.formula(fs), data=eval(mcall$data))[,1]
        spaceref <- 1
        spaceref <- unique(eval(mcall$control.st$spaceref))
        if (length(spaceref)>1)
            stop("length(spaceref)>1")
        if (length(spaceref)==0) spaceref <- 1
    }
    fst <- paste('~', mcall$control.st$spacetime)
    if (fst=='~ ') {
        spacetime <- NULL
    } else {
        stname <- substring(fst, 3)
        spacetime <- model.frame(as.formula(fst), data=eval(mcall$data))[,1]
    }
###    cat('tname =', tname, ', sname =', sname, ', stname =', stname, '\n')
    type <- as.character(unique(eval(mcall$control.st$type)))
    types <- c(1:4, paste0(2:4, 'c'), paste0(2:4, 'd'))
    type <- if(length(type)==0) types else types[match(type, types)]
###    cat('type =', type, '\n')
    if (length(type)==0) return(NULL)
    else res <- list()
    if (length(mcall$control.st$diagonal)==0) diagonal <- 1e-5
###    cat('diagonal =', diagonal, '\n')
    m <- n <- NULL
    nst <- length(unique(spacetime))
    if (!is.null(mcall$control.st$graph)) {
        graph <- eval(mcall$control.st$graph)
        n <- nrow(graph <- inla.graph2matrix(graph))
        R.s <- inla.scale.model(Diagonal(n, colSums(graph)) - graph,
                                constr=list(A=matrix(1, 1, n), e=0))
    } else {
        if (any(substr(type,1,1)%in%c('3', '4')))
            stop("'graph' must be provided to build the spacetime interaction model!")
    }
    if (!is.null(space)) {
        if (!is.null(mcall$control.st$graph))
            if(n!=length(unique(space)))
                stop("Size of 'space' is not equal to the size of 'graph'!")
        n <- length(unique(space))
    }
    if (!is.null(time)) {
        m <- length(unique(time))
        if (any(substr(type,1,1)%in%c('2', '4')))
            R.t <- inla.scale.model(crossprod(diff(Diagonal(m))),
                                   constr=list(A=matrix(1,1,m), e=0))
    }
    if (is.null(m)) m <- nst/n
    if (is.null(n)) n <- nst/m
###    cat('m = ', m, ', n = ', n, ', nst = ', nst, '\n', sep='')
    if (TRUE) { ## working in progress: identify need of constraints from the formula
        etemp <- INLA:::inla.interpret.formula(formula, data, debug=FALSE)
        rterms <- attr(terms(etemp[[1]]), 'term.labels')
        id.r <- which(sapply(etemp$random.spec, function(x) is.null(x$weights))) 
        if (length(id.r)>0) { 
            r.rankdef <- which(sapply(etemp$random.spec[id.r], function(x) is.null(x$rankdef))) 
            if (length(r.rankdef)>0) { 
                r.size <- sapply(etemp$random.spec[id.r[r.rankdef]], function(x) x$n) 
                j.s <- which(r.size==n) 
                j.t <- which(r.size==m) 
                if (length(j.s)>1) 
                    stop('Too many spatial effects with rank deficiency.') 
                if (length(j.t)>1) 
                    stop('Too many temporal effects with rank deficiency.') 
            }
        }
        lc2.on <- any(rterms==sname)
        lc3.on <- any(rterms==tname)
    }
###    cat('lc2 =', lc2.on, ' and lc3 =', lc3.on, '\n')
    M2 <- kronecker(matrix(1/m,1,m), diag(n)) 
    M3 <- kronecker(diag(m), matrix(1/n,1,n))
    dotdot <- mcall$control.st[which(!is.element(names(mcall$control.st),
                                                 c('time', 'space', 'graph', 'type',
                                                   'timeref', 'spaceref')))]

    add0 <- ''
    if(length(names(dotdot))>2) {
        add0 <- paste0(', ', names(dotdot)[3], '=', dotdot[3])
    }
    if (any(type%in%'1')) {
        add1 <- paste0('f(', stname, ', model="iid"', add0, ')')
        res$'1' <- inla(update(formula, paste('.~.+', add1)), ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='1') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'2')) {
        add2 <- paste0('f(', stname, ', model="generic0", constr=FALSE, ',
                       'Cmatrix=kronecker(R.t, Diagonal(n)), ',
                       'extraconstr=list(A=M2, e=rep(0,n))', add0, ')')
        res$'2' <- inla(update(formula, paste('.~.+', add2)), ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='2') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'3')) {
        add3 <- paste0('f(', stname, ', model="generic0", constr=FALSE, ',
                       'Cmatrix=kronecker(Diagonal(m), R.s), ',
                       'extraconstr=list(A=M3, e=rep(0,m))', add0, ')')
        res$'3' <- inla(update(formula, paste('.~.+', add3)), ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='3') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'4')) {
        add4 <- paste0('f(', stname, ', model="generic0", constr=FALSE, ',
                       'Cmatrix=kronecker(R.t, R.s), ',
                       'extraconstr=list(A=rbind(M2,M3), e=rep(0,n+m))', add0, ')')
        res$'4' <- inla(update(formula, paste('.~.+', add4)), ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='4') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'2c')) {
        add2c <- paste0('f(st2, model="generic0", constr=FALSE, ',
                        'Cmatrix=kronecker(R.t, Diagonal(n)), ',
                        'extraconstr=list(A=M2, e=rep(0,n))', add0, ')')
        if(is.null(time)) {
            st2 <- spacetime
        } else {
            st2 <- ifelse(time==timeref, NA, spacetime-cumsum(time==timeref)) 
        }
        id2 <- which(!is.na(st2))
        lc2 <- inla.make.lincombs(
            st2=Diagonal(n*m)[,id2] - M2[space, id2])
        names(lc2) <- gsub('lc', 'st', names(lc2))
        if (lc2.on) {
            lc2args <- list(cbind(Diagonal(n), Diagonal(n,0)), st2=M2[,id2]-1/(n*m))
            names(lc2args)[1] <- sname 
            lcc2 <- do.call('inla.make.lincombs', lc2args) 
            names(lcc2) <- gsub('lc', 's', names(lcc2))
            lc2 <- c(lcc2, lc2)
        } 
        res$'2c' <- inla(update(formula, paste('.~.+', add2c)),
                         lincomb=lc2, ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='2c') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'3c')) {
        add3c <- paste0('f(st3, model="generic0", constr=FALSE, ',
                        'Cmatrix=kronecker(Diagonal(m), R.s), ',
                        'extraconstr=list(A=M3, e=rep(0,m))', add0, ')')
        if(is.null(space)) {
            st3 <- spacetime
        } else {
            st3 <- ifelse(space==spaceref, NA, spacetime-cumsum(space==spaceref)) 
        }
        id3 <- which(!is.na(st3))
        lc3 <- inla.make.lincombs(
            st3=Diagonal(n*m)[,id3] - M3[time, id3])
        names(lc3) <- gsub('lc', 'st', names(lc3))
        if (lc3.on) {
            lc3args <- list(Diagonal(m), st3=M3[,id3]-1/(n*m))
            names(lc3args)[1] <- tname
            lcc3 <- do.call('inla.make.lincombs', lc3args) 
            names(lcc3) <- gsub('lc', 't', names(lcc3))
            lc3 <- c(lcc3, lc3)
        }
        res$'3c' <- inla(update(formula, paste('.~.+', add3c)),
                         lincomb=lc3, ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='3c') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'4c')) {
        add4c <- paste0('f(st4, model="generic0", constr=FALSE, ',
                        'Cmatrix=kronecker(R.t, R.s), ',
                        'extraconstr=list(A=rbind(M2, M3), e=rep(0,n+m))', add0, ')')
        aux4 <- logical(n*m) 
        if (is.null(time)) 
            aux4 <- aux4|(time==timeref)
        if (is.null(space))
            aux4 <- aux4|(space==spaceref)
        st4 <- ifelse(aux4, NA, spacetime-cumsum(aux4)) 
        id4 <- which(!is.na(st4))
        lc4 <- inla.make.lincombs(
            st4=(Diagonal(n*m)[,id4] - M2[space,id4] -M3[time,id4] +1/(n*m)))
        names(lc4) <- gsub('lc', 'st', names(lc4))
        if (lc3.on) {
            lc3args <- list(Diagonal(m), st4=M3[,id4]-1/(n*m))
            names(lc3args)[1] <- tname
            lcc3 <- do.call('inla.make.lincombs', lc3args)
            names(lcc3) <- gsub('lc', 't', names(lcc3))
            lc4 <- c(lcc3, lc4)
        }
        if (lc2.on) {
            lc2args <- list(cbind(Diagonal(n), Diagonal(n,0)),
                            st4=M2[, id4]-1/(n*m))
            names(lc2args)[1] <- sname
            lcc2 <- do.call('inla.make.lincombs', lc2args)
            names(lcc2) <- gsub('lc', 's', names(lcc2))
            lc4 <- c(lcc2, lc4)
        }
        res$'4c' <- inla(update(formula, paste('.~.+', add4c)),
                         lincomb=lc4, ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='4c') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%c('2d', '3d', '4d'))) {
        lcd2 <- lcd3 <- NULL
        if (lc2.on) {
            lc2args <- list(cbind(Diagonal(n), Diagonal(n,0)), M2)
            names(lc2args) <- c(sname, stname)
            lcd2 <- do.call('inla.make.lincombs', lc2args)
        }
        if (lc3.on) {
            lc3args <- list(Diagonal(m), M3)
            names(lc3args) <- c(tname, stname)
            lcd3 <- do.call('inla.make.lincombs', lc3args)
        }
        dd <- Diagonal(m*n, diagonal) 
    }
    if (any(type%in%'2d')) {
        lcd2args <- list(Diagonal(n*m)-M2[space,])
        names(lcd2args) <- stname
        lcd <- do.call('inla.make.lincombs', lcd2args) 
        add2d <- paste0('f(', stname, ', model="generic0", ',
                        'constr=TRUE, rankdef=n, ',
                        'Cmatrix=kronecker(R.t, Diagonal(n)) + dd', add0, ')')
        res$'2d' <- inla(update(formula, paste('.~.+', add2d)),
                         lincomb=c(lcd2, lcd), ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='2d') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'3d')) {
        lcd3args <- list(Diagonal(n*m) -M3[time,])
        names(lcd3args) <- stname
        lcd <- do.call('inla.make.lincombs', lcd3args) 
        names(lcd) <- gsub('lc', 'st', names(lcd))
        add3d <- paste0('f(', stname, ', model="generic0", ',
                        'constr=TRUE, rankdef=m, ',
                        'Cmatrix=kronecker(Diagonal(m), R.s) + dd', add0, ')')
        res$'3d' <- inla(update(formula, paste('.~.+', add3d)), 
                         lincomb=c(lcd3, lcd), ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='3d') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (any(type%in%'4d')) {
        lcd3args <- lcd2args <- NULL
        if(lc2.on) {
            lcd2args <- list(cbind(Diagonal(n), Diagonal(n,0)), M2)
            names(lcd2args) <- c(sname, stname)
            lcd2 <- do.call('inla.make.lincombs', lcd2args)
            names(lcd2) <- gsub('lc', 's', names(lcd2))
        }
        if(lc3.on) {
            lcd3args <- list(Diagonal(m), M3)
            names(lcd3args) <- c(tname, stname)
            lcd3 <- do.call('inla.make.lincombs', lcd3args)
            names(lcd3) <- gsub('lc', 't', names(lcd3))
        }
        lcdargs <- list(Diagonal(n*m) - M2[space,] - M3[time,] +1/(n*m))
        names(lcdargs) <- stname
        lcd <- do.call('inla.make.lincombs', lcdargs)
        names(lcd) <- gsub('lc', 'st', names(lcd))
        add4d <- paste0('f(', stname, ', model="generic0", ', 
                        'constr=TRUE, rankdef=n+m, ', 
                        'Cmatrix=kronecker(R.t, R.s) + dd', add0, ')')
        res$'4d' <- inla(update(formula, paste('.~.+', add4d)),
                         lincomb=c(lcd2, lcd3, lcd), ...)
    }
    if(progress && (length(res)>0) && tail(names(res),1)=='4d') 
        cat('type = ', tail(names(res),1), ', cpu = ',  res[[length(res)]]$cpu[4], '\n', sep='')
    if (length(res)==1) return(res[[1]]) 
    return(res) 
}


