#' Spacetime interaction models
#' 
#' It implements the models in Knorr-Held, L. (2000) with three different
#' constraint approaches: sum-to-zero, contrast or diagonal add.
#' 
#' @rdname knmodels
#' @aliases inla.knmodels knmodels
#' @param formula The formula specifying the other model components, without
#' the spacetime interaction term. The spacetime interaction term will be added
#' accordly to the specification in the `control.st` argument. See
#' `inla`
#' @param progress If it is to be shown the model fitting progress. Useful if
#' more than one interaction type is being fitted.
#' @param control.st Named list of arguments to control the spacetime
#' interaction. It should contain:
#' \describe{
#'   \item{time}{to be used as the
#' index set for the main temporal effect which will be considered for the
#' constraints when it is the case.}
#'   \item{space}{to be used as the index set
#' for the main spatial effect which will be considered for the constraints
#' when it is the case.}
#'   \item{spacetime}{to be the index set for the spacetime
#' interaction effect.}
#'   \item{graph}{to be the graph for the spatial neighbor
#' structure to be used in a [f()] term for the main spatial random
#' effect term or for building the spacetime interaction model.}
#'   \item{type}{to
#' specify the spacetime interaction type.  `1` to `4` corresponds to
#' the four interaction types in Knorr-Held, L. (2000) with all the needed
#' sum-to-zero constraints.  `2c`, `3c` and `4c` are the
#' contrast version considering the first time or space constrained to be equal
#' to zero.  `2d`, `3d` and `4d` are the corresponding versions
#' when considering the diagonal add approach.}
#'   \item{diagonal}{to be the value
#' to be added to the diagonal when using the diagonal add approach.}
#'   \item{timeref}{to specify the time point to be the reference time in the
#' contrast parametrization.}
#'   \item{spaceref}{to specify the area to be the
#' reference for the contrast parametrization.}
#'   \item{...}{where additional
#' arguments can be passed to [f()] function.  Specification of the
#' hyperparameter, fixed or random, initial value, prior and its parameters for
#' the spacetime interaction. See `?inla.models` and look for
#' `generic0`.  By default we scale it and use the PC-prior to set the
#' prior using the `pc.prec` prior with `param = c(0.5, 0.5)`. See
#' documentation with `?inla.doc("pc.prec")`.} }
#' @param ... Arguments to be passed to the [inla()] function.
#' @param envir Environment in which to evaluate the ... arguments.
#' @return `inla.knmodels` returns an object of class `"inla"`.  or a
#' list of objects of this class if it is asked to compute more than one
#' interaction type at once.  Note: when the model type is 2c, 3c, 4c, 2d, 3d
#' or 4d, it also includes linear combinations summary.
#' @author Elias T. Krainski
#' @seealso [inla.knmodels.sample()] to sample from
#' @examplesIf require("sp")
#' 
#' ### define space domain as a grid
#' grid <- sp::SpatialGrid(sp::GridTopology(c(0,0), c(1, 1), c(4, 5)))
#' (n <- nrow(xy <- sp::coordinates(grid)))
#' 
#' ### build a spatial neighborhood list
#' jj <- lapply(1:n, function(i)
#'     which(sqrt((xy[i,1]-xy[,1])^2 + (xy[i,2]-xy[,2])^2)==1))
#' 
#' ### build the spatial adjacency matrix
#' graph <- sparseMatrix(rep(1:n, sapply(jj, length)),
#'                       unlist(jj), x=1, dims=c(n, n))
#' 
#' ### some random data at 10 time point
#' dat <- inla.knmodels.sample(graph, m=10, tau.t=2, tau.s=2, tau.st=3)
#' str(dat)
#' sapply(dat$x, summary)
#' 
#' nd <- length(dat$x$eta)
#' dat$e <- runif(nd, 0.9, 1.1)*rgamma(n, 40, 2)
#' dat$y <- rpois(nd, dat$e*exp(dat$x$eta-3))
#' summary(dat$y)
#'
#' ### fit the type 4 considering three different approaches
#' tgraph <- sparseMatrix(i=c(2:10, 1:9), j=c(1:9, 2:10), x=1)
#' res <- inla.knmodels(y ~ f(time, model='bym2', graph=tgraph) +
#'      f(space, model='bym2', graph=graph),
#'      data=dat, family='poisson', E=dat$E, progress=TRUE,
#'      control.st=list(time=time, space=space,
#'         spacetime=spacetime, graph=graph, type=c(4, '4c')), 
#'      control.compute=list(dic=TRUE, waic=TRUE, cpo=TRUE))
#' sapply(res, function(x)
#'        c(dic=x$dic$dic, waic=x$waic$waic, cpo=-sum(log(x$cpo$cpo))))
#' 
#' @export inla.knmodels
`inla.knmodels` <- function(formula,
                            progress = FALSE,
                            control.st = list(
                                time,
                                space,
                                spacetime,
                                graph,
                                type = c(paste(1:4), paste0(2:4, "c"), paste0(2:4, "d")),
                                diagonal = 1e-5,
                                timeref = 1,
                                spaceref = 1
                            ),
                            ...,
                            envir = parent.frame()) {
    mcall <- match.call(expand.dots = TRUE)
    ft <- paste("~", mcall$control.st$time)
    if (ft == "~ ") {
        time <- tname <- NULL
    } else {
        tname <- substring(ft, 3)
        time <- model.frame(as.formula(ft), data = eval(mcall$data, envir = envir))[, 1]
        timeref <- unique(eval(mcall$control.st$timeref, envir = envir))
        if (length(timeref) > 1) {
            stop("length(timeref)>1")
        }
        if (length(timeref) == 0) timeref <- 1
    }
    fs <- paste("~", mcall$control.st$space)
    if (fs == "~ ") {
        space <- sname <- NULL
    } else {
        sname <- substring(fs, 3)
        space <- model.frame(as.formula(fs), data = eval(mcall$data, envir = envir))[, 1]
        spaceref <- 1
        spaceref <- unique(eval(mcall$control.st$spaceref, envir = envir))
        if (length(spaceref) > 1) {
            stop("length(spaceref)>1")
        }
        if (length(spaceref) == 0) spaceref <- 1
    }
    assign('timeref', timeref,
           envir=environment(formula))
    assign('spaceref', spaceref,
           envir=environment(formula))
    fst <- paste("~", mcall$control.st$spacetime)
    if (fst == "~ ") {
        spacetime <- NULL
    } else {
        stname <- substring(fst, 3)
        spacetime <- model.frame(as.formula(fst), data = eval(mcall$data, envir = envir))[, 1]
    }
    ## cat('tname =', tname, ', sname =', sname, ', stname =', stname, '\n')
    type <- as.character(unique(eval(mcall$control.st$type, envir = envir)))
    types <- c(1:4, paste0(2:4, "c"), paste0(2:4, "d"))
    type <- if (length(type) == 0) types else types[match(type, types)]
    ## cat('type =', type, '\n')
    if (length(type) == 0) {
        return(NULL)
    } else {
        res <- list()
    }
    if (length(mcall$control.st$diagonal) == 0) diagonal <- 1e-5
    ## cat('diagonal =', diagonal, '\n')
    .no.of.t <- .no.of.s <- NULL
    nst <- length(unique(spacetime))
    if (!is.null(mcall$control.st$graph)) {
        graph <- eval(mcall$control.st$graph, envir = envir)
        .no.of.s <- nrow(graph <- inla.graph2matrix(graph))
        R.s <- inla.scale.model(Diagonal(.no.of.s, colSums(graph)) - graph,
            constr = list(A = matrix(1, 1, .no.of.s), e = 0)
        )
        assign('R.s', R.s,
               envir=environment(formula))
    } else {
        if (any(substr(type, 1, 1) %in% c("3", "4"))) {
            stop("'graph' must be provided to build the spacetime interaction model!")
        }
    }
    if (!is.null(space)) {
        if (!is.null(mcall$control.st$graph)) {
            if (nrow(graph) != length(unique(space))) {
                stop("Size of 'space' is not equal to the size of 'graph'!")
            }
        }
        n <- length(unique(space))
    }
    if (!is.null(time)) {
        .no.of.t <- length(unique(time))
        if (any(substr(type, 1, 1) %in% c("2", "4"))) {
            R.t <- inla.scale.model(
                       crossprod(diff(Diagonal(.no.of.t))),
                constr = list(A = matrix(1, 1, .no.of.t), e = 0))
            assign('R.t', R.t,
                   envir=environment(formula))
        }
    }
    if (is.null(.no.of.t)) .no.of.t <- nst / .no.of.s
    if (is.null(.no.of.s)) .no.of.s <- nst / .no.of.t

    assign('.no.of.s', .no.of.s,
           envir=environment(formula))
    assign('.no.of.t', .no.of.t,
           envir=environment(formula))

##     cat('.no.of.t = ', .no.of.t, ', .no.of.s = ', .no.of.s, ', nst = ', nst, '\n', sep='')
    if (TRUE) { ## working in progress: identify need of constraints from the formula
        etemp <- inla.interpret.formula(formula, data, debug = FALSE)
##        print(str(etemp))
        rterms <- attr(terms(etemp[[1]]), "term.labels")
        if(length(rterms)>0) {
            r.size <- sapply(rterms, function(x)
                length(unique(eval(mcall$data, envir = envir)[[x]])))
  ##          print(r.size)
            r.size.r <- sapply(etemp$random.spec, function(x)
            (!is.null(x$rankdef)) |
            (x$model %in% c("rw1", "rw2", "besag", 'bym', 'bym2'))) * r.size
    ##          print(r.size.r)
            j.s <- which(r.size.r == .no.of.s)
            j.t <- which(r.size.r == .no.of.t)
            if (length(j.s) > 1) {
                stop("Too many spatial effects with rank deficiency.")
            }
            if (length(j.t) > 1) {
                stop("Too many temporal effects with rank deficiency.")
            }
        }
        lc2.on <- any(rterms == sname)
        lc3.on <- any(rterms == tname)
	if(progress) {
		cat('r.size.r =', r.size.r, '\n')
	}
    }
    ## cat('lc2 =', lc2.on, ' and lc3 =', lc3.on, '\n')
    M2 <- kronecker(matrix(1 / .no.of.t, 1, .no.of.t), diag(.no.of.s))
    M3 <- kronecker(diag(.no.of.t), matrix(1 / .no.of.s, 1, .no.of.s))
    assign('M2', M2,
           envir=environment(formula))
    assign('M3', M3,
           envir=environment(formula))

    dotdot <- mcall$control.st[which(!is.element(
        names(mcall$control.st),
        c(
            "time", "space", "graph", "type",
            "timeref", "spaceref"
        )
    ))]

    add0 <- ""
    if (length(names(dotdot)) > 2) {
        add0 <- paste0(", ", names(dotdot)[3], "=", dotdot[3])
    }
    if (any(type %in% "1")) {
        add1 <- paste0("f(", stname, ', model="iid"', add0, ")")
        res$"1" <- inla(update(formula, paste(".~.+", add1)), ...)
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "1") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "2")) {
        add2 <- paste0(
            "f(", stname, ', model="generic0", constr=FALSE, ',
            "Cmatrix=kronecker(R.t, Diagonal(.no.of.s))",
            ifelse(lc2.on, ", extraconstr=list(A=M2, e=rep(0,.no.of.s))", ""), 
	    add0, ")"
        )
	if(progress) 
	  cat("Added model term:\n", gsub(".no.of.s", .no.of.s, add2, fixed = TRUE), "\n")
        res$"2" <- inla(update(formula, paste(".~.+", add2)), ...)
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "2") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "3")) {
        add3 <- paste0(
            "f(", stname, ', model="generic0", constr=FALSE, ',
            "Cmatrix=kronecker(Diagonal(.no.of.t), R.s)",
            ifelse(lc3.on, ", extraconstr=list(A=M3, e=rep(0,.no.of.t))", ""), 
	    add0, ")"
        )
        if(progress)
          cat("Added model term:\n", gsub(".no.of.t", .no.of.t, add3, fixed = TRUE), "\n")
        res$"3" <- inla(update(formula, paste(".~.+", add3)), ...)
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "3") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "4")) {
	if(lc2.on) {
	  if(lc3.on) {
            add4 <- paste0(
             "f(", stname, ', model="generic0", constr=FALSE, ',
             "Cmatrix=kronecker(R.t, R.s), ",
             "extraconstr=list(A=rbind(M2,M3)[-1,], e=rep(0,.no.of.s+.no.of.t-1))", add0, ")"
	     )
	  } else {
            add4 <- paste0(
             "f(", stname, ', model="generic0", constr=FALSE, ',
             "Cmatrix=kronecker(R.t, R.s), ",
             "extraconstr=list(A=M2, e=rep(0,.no.of.s))", add0, ")"
           )
	  }
        } else {
	  if(lc3.on) {
           add4 <- paste0(
             "f(", stname, ', model="generic0", constr=FALSE, ',
             "Cmatrix=kronecker(R.t, R.s), ",
             "extraconstr=list(A=M3, e=rep(0,.no.of.t))", add0, ")"
             )
          } else {
            add4 <- paste0(
             "f(", stname, ', model="generic0", constr=FALSE, ',
             "Cmatrix=kronecker(R.t, R.s)", add0, ")"
           )
          }
	}
        if(progress) {
          cat("Added model term:\n", 
	      gsub(".no.of.s", .no.of.s, 
		   gsub(".no.of.t", .no.of.t, add4, fixed = TRUE), 
		   fixed = TRUE), "\n")
	}
        res$"4" <- inla(update(formula, paste(".~.+", add4)), ...)
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "4") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "2c")) {
        add2c <- paste0(
            'f(st2, model="generic0", constr=FALSE, ',
            "Cmatrix=kronecker(R.t[-timeref, -timeref], Diagonal(.no.of.s))", add0, ")"
        )
        if (is.null(time)) {
            st2 <- spacetime
        } else {
            st2 <- ifelse(time == timeref, NA, spacetime - cumsum(time == timeref))
        }
        assign('st2', st2, environment(formula))
        id2 <- which(!is.na(st2))
        lc2 <- inla.make.lincombs(
            st2 = Diagonal(.no.of.s * .no.of.t)[, id2] - M2[space, id2]
        )
        names(lc2) <- gsub("lc", "st", names(lc2))
        if (lc2.on) {
            lc2args <- list(cbind(Diagonal(.no.of.s), Diagonal(.no.of.s, 0)), st2 = M2[, id2] - 1 / (.no.of.s * .no.of.t))
            names(lc2args)[1] <- sname
            lcc2 <- do.call("inla.make.lincombs", lc2args)
            names(lcc2) <- gsub("lc", "s", names(lcc2))
            lc2 <- c(lcc2, lc2)
        }
        res$"2c" <- inla(update(formula, paste(".~.+", add2c)),
            lincomb = lc2, ...
        )
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "2c") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "3c")) {
        add3c <- paste0(
            'f(st3, model="generic0", constr=FALSE, ',
            "Cmatrix=kronecker(Diagonal(.no.of.t), R.s[-spaceref, -spaceref])", add0, ")"
        )
        if (is.null(space)) {
            st3 <- spacetime
        } else {
            st3 <- ifelse(space == spaceref, NA, spacetime - cumsum(space == spaceref))
        }
        assign('st3', st3, environment(formula))
        id3 <- which(!is.na(st3))
        lc3 <- inla.make.lincombs(
            st3 = Diagonal(.no.of.s * .no.of.t)[, id3] - M3[time, id3]
        )
        names(lc3) <- gsub("lc", "st", names(lc3))
        if (lc3.on) {
            lc3args <- list(Diagonal(.no.of.t), st3 = M3[, id3] - 1 / (.no.of.s * .no.of.t))
            names(lc3args)[1] <- tname
            lcc3 <- do.call("inla.make.lincombs", lc3args)
            names(lcc3) <- gsub("lc", "t", names(lcc3))
            lc3 <- c(lcc3, lc3)
        }
        res$"3c" <- inla(update(formula, paste(".~.+", add3c)),
            lincomb = lc3, ...
        )
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "3c") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "4c")) {
        add4c <- paste0(
            'f(st4, model="generic0", constr=FALSE, ',
            "Cmatrix=kronecker(R.t[-timeref, -timeref], R.s[-spaceref, -spaceref])", add0, ")"
        )
        if (is.null(time)) {
            if(is.null(space)) {
                st4 <- spacetime
            } else {
                st4 <- ifelse(space == spaceref, NA, spacetime - cumsum(space == spaceref))
            }
        } else {
            if(is.null(space)) {
                st4 <- ifelse(time == timeref, NA, spacetime - cumsum(time == timeref))
            } else {
                st4 <- ifelse(
                    time == timeref, NA,
                       ifelse(space == spaceref, NA,
                              spacetime - cumsum((time == timeref) |
                                                 (space == spaceref))))
            }
        }
        assign('st4', st4, environment(formula))
        id4 <- which(!is.na(st4))
        lc4 <- inla.make.lincombs(
            st4 = (Diagonal(.no.of.s * .no.of.t)[, id4] - M2[space, id4] - M3[time, id4] + 1 / (.no.of.s * .no.of.t))
        )
        names(lc4) <- gsub("lc", "st", names(lc4))
        if (lc3.on) {
            lc3args <- list(Diagonal(.no.of.t), st4 = M3[, id4] - 1 / (.no.of.s * .no.of.t))
            names(lc3args)[1] <- tname
            lcc3 <- do.call("inla.make.lincombs", lc3args)
            names(lcc3) <- gsub("lc", "t", names(lcc3))
            lc4 <- c(lcc3, lc4)
        }
        if (lc2.on) {
            lc2args <- list(cbind(Diagonal(.no.of.s), Diagonal(.no.of.s, 0)),
                st4 = M2[, id4] - 1 / (.no.of.s * .no.of.t)
            )
            names(lc2args)[1] <- sname
            lcc2 <- do.call("inla.make.lincombs", lc2args)
            names(lcc2) <- gsub("lc", "s", names(lcc2))
            lc4 <- c(lcc2, lc4)
        }
        res$"4c" <- inla(update(formula, paste(".~.+", add4c)),
            lincomb = lc4, ...
        )
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "4c") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% c("2d", "3d", "4d"))) {
        lcd2 <- lcd3 <- NULL
        if (lc2.on) {
            lc2args <- list(cbind(Diagonal(.no.of.s), Diagonal(.no.of.s, 0)), M2)
            names(lc2args) <- c(sname, stname)
            lcd2 <- do.call("inla.make.lincombs", lc2args)
        }
        if (lc3.on) {
            lc3args <- list(Diagonal(.no.of.t), M3)
            names(lc3args) <- c(tname, stname)
            lcd3 <- do.call("inla.make.lincombs", lc3args)
        }
        dd <- Diagonal(.no.of.t * .no.of.s, diagonal)
        assign('dd', dd, environment(formula))
    }
    if (any(type %in% "2d")) {
        lcd2args <- list(Diagonal(.no.of.s * .no.of.t) - M2[space, ])
        names(lcd2args) <- stname
        lcd <- do.call("inla.make.lincombs", lcd2args)
        add2d <- paste0(
            "f(", stname, ', model="generic0", ',
            "constr=TRUE, rankdef=.no.of.s, ",
            "Cmatrix=kronecker(R.t, Diagonal(.no.of.s)) + dd", add0, ")"
        )
        res$"2d" <- inla(update(formula, paste(".~.+", add2d)),
            lincomb = c(lcd2, lcd), ...
        )
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "2d") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "3d")) {
        lcd3args <- list(Diagonal(.no.of.s * .no.of.t) - M3[time, ])
        names(lcd3args) <- stname
        lcd <- do.call("inla.make.lincombs", lcd3args)
        names(lcd) <- gsub("lc", "st", names(lcd))
        add3d <- paste0(
            "f(", stname, ', model="generic0", ',
            "constr=TRUE, rankdef=.no.of.t, ",
            "Cmatrix=kronecker(Diagonal(.no.of.t), R.s) + dd", add0, ")"
        )
        res$"3d" <- inla(update(formula, paste(".~.+", add3d)),
            lincomb = c(lcd3, lcd), ...
        )
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "3d") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (any(type %in% "4d")) {
        lcd3args <- lcd2args <- NULL
        if (lc2.on) {
            lcd2args <- list(cbind(Diagonal(.no.of.s), Diagonal(.no.of.s, 0)), M2)
            names(lcd2args) <- c(sname, stname)
            lcd2 <- do.call("inla.make.lincombs", lcd2args)
            names(lcd2) <- gsub("lc", "s", names(lcd2))
        }
        if (lc3.on) {
            lcd3args <- list(Diagonal(.no.of.t), M3)
            names(lcd3args) <- c(tname, stname)
            lcd3 <- do.call("inla.make.lincombs", lcd3args)
            names(lcd3) <- gsub("lc", "t", names(lcd3))
        }
        lcdargs <- list(Diagonal(.no.of.s * .no.of.t) - M2[space, ] - M3[time, ] + 1 / (.no.of.s * .no.of.t))
        names(lcdargs) <- stname
        lcd <- do.call("inla.make.lincombs", lcdargs)
        names(lcd) <- gsub("lc", "st", names(lcd))
        add4d <- paste0(
            "f(", stname, ', model="generic0", ',
            "constr=TRUE, rankdef=.no.of.s+.no.of.t, ",
            "Cmatrix=kronecker(R.t, R.s) + dd", add0, ")"
        )
        res$"4d" <- inla(update(formula, paste(".~.+", add4d)),
            lincomb = c(lcd2, lcd3, lcd), ...
        )
    }
    if (progress && (length(res) > 0) && tail(names(res), 1) == "4d") {
        cat("type = ", tail(names(res), 1), ", cpu = ", res[[length(res)]]$cpu.used[4], "\n", sep = "")
    }
    if (length(res) == 1) {
        return(res[[1]])
    }
    return(res)
}
