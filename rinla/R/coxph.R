#' Convert a Cox proportional hazard model into Poisson regression
#' 
#' Tools to convert a Cox proportional hazard model into Poisson regression
#' 
#' 
#' @aliases inla.coxph inla.rbind.data.frames coxph cbind.data.frames
#' @param formula The formula for the coxph model where the response must be a
#' `inla.surv`-object.
#' @param data All the data used in the formula, as a list.
#' @param control.hazard Control the model for the baseline-hazard; see
#' `?control.hazard`.
#' @param tag An optional tag added to the names of the new variables created
#' (to make them unique when combined with several calls of `inla.coxph`.
#' Note that `E..coxph` is not included, as its usually merged into one
#' vector over different expansions.
#' @param debug Print debug-information
#' @param ... Data.frames to be `rbind`-ed, padding with `NA`.
#' @return `inla.coxph` returns a list of new expanded variables to be
#' used in the `inla`-call.  Note that element `data` and
#' `data.list` needs to be merged into a `list` to be passed as the
#' `data` argument. See the example for details.
#' 
#' `inla.rbind.data.frames` returns the rbinded data.frames padded with
#' NAs.  There is a better implementation in `dplyr::bind_rows`, which is
#' used if package `dplyr` is installed.
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @examples
#' 
#' ## How the cbind.data.frames works:
#' df1 = data.frame(x=1:2, y=2:3, z=3:4)
#' df2 = data.frame(x=3:4, yy=4:5, zz=5:6)
#' inla.rbind.data.frames(df1, df2)
#' 
#' ## Standard example of how to convert a coxph into a Poisson regression
#' n = 1000
#' x = runif(n)
#' lambda = exp(1+x)
#' y = rexp(n, rate=lambda)
#' event = rep(1,n)
#' data = list(y=y, event=event, x=x)
#' y.surv = inla.surv(y, event)
#' intercept1 = rep(1, n)
#' p = inla.coxph(y.surv ~ -1 + intercept1 + x,
#'                list(y.surv = y.surv,  x=x, intercept1 = intercept1))
#' 
#' r = inla(p$formula,
#'         family = p$family,
#'         data=c(as.list(p$data), p$data.list),
#'         E = p$E)
#' summary(r)
#' 
#' ## How to use this in a joint model
#' intercept2 = rep(1, n)
#' y = 1 + x + rnorm(n, sd=0.1)
#' df = data.frame(intercept2, x, y)
#' 
#' ## new need to cbind the data.frames, and then add the list-part of
#' ## the data
#' df.joint = c(as.list(inla.rbind.data.frames(p$data, df)), p$data.list)
#' df.joint$Y = cbind(df.joint$y..coxph, df.joint$y)
#' 
#' ## merge the formulas, recall to add '-1' and to use the new joint
#' ## reponse 'Y'
#' formula = update(p$formula, Y ~ intercept2 -1 + .)
#' 
#' rr = inla(formula,
#'         family = c(p$family, "gaussian"),
#'         data = df.joint,
#'         E = df.joint$E..coxph)
#' 
#' @rdname coxph
#' @export inla.coxph
`inla.coxph` <- function(formula, data, control.hazard = list(), debug = FALSE, tag = "") 
{
    ## convert a coxph-model into poisson-regression and return a new
    ## data-list with the expand variables and new variables to use in
    ## the poisson regression

    if (is.data.frame(data)) {
        data <- as.list(data)
    }

    name.y <- inla.formula2character(formula[2])
    tmp <- (names(data) %in% name.y)
    y.surv <- NULL
    if (any(tmp)) {
        if (sum(tmp) > 1) {
            stop(paste(c("Several entries in 'data' match the name of the response:",
                         "response=", name.y, ", matches=", sum(tmp), ".")))
        }
        y.surv <- data[[which(tmp)]]
        data[[which(tmp)]] <- NULL
    } else {
        try.res <- try(inla.eval(paste(c("y.surv = with(data,", name.y, ")"))), silent = TRUE)
        if (inherits(try.res, "try-error")) {
            stop(paste(c("The reponse '", name.y, "' is not in 'data' and trying to expand it, failed.")))
        }
    }

    len.y.surv <- max(sapply(y.surv, length))
    data.f <- inla.fix.data(data, len.y.surv, revert = FALSE)
    data.l <- inla.fix.data(data, len.y.surv, revert = TRUE)
    data.f <- try(as.data.frame(data.f), silent = TRUE)
    if (inherits(data.f, "try-error")) {
        stop("Failed to convert 'data' into a 'data.frame' even after 'inla.fix.data'.")
    }

    if (!inherits(y.surv, "inla.surv")) {
        stop(paste0("For survival models, then the reponse has to be of class `inla.surv'; you have `",
                    paste0(class(y.surv), collapse = ", "),
                    "'"))
    }


    control.hazard <- ctrl_object(control.hazard, "hazard", data.f)
    cont.hazard <- ctrl_update(control.hazard)

    if (is.null(y.surv$subject)) {
        res <- inla.expand.dataframe.1(y.surv, data.f, control.hazard = cont.hazard)
    } else {
        res <- inla.expand.dataframe.2(y.surv, data.f, control.hazard = cont.hazard)
    }

    strata.var <- NULL
    strata.tmp <- NULL
    if (!is.null(cont.hazard$strata.name)) {
        if (is.character(cont.hazard$strata.name) && length(cont.hazard$strata.name) == 1) {
            ## strata = "x"
            strata.var <- cont.hazard$strata
        } else {
            stop("Argument to `strata.name' must be the name of a variable in the data.frame.")
        }
    }
    if (!is.null(strata.var)) {
        if (!is.element(strata.var, names(data.f))) {
            stop(inla.paste(c("Variable `", strata.var,
                              "' in control.hazard=list(strata=...) needs to be in the data.frame: names(data) = ",
                              names(data.f))))
        }
        if (debug) print("apply inla.strata() on strata.var")
        ## cleaner code...
        strata.tmp <- inla.strata(inla.get.var(strata.var, data.f))
        data.f[[strata.var]] <- strata.tmp$strata
        data.l$baseline.hazard.strata.coding <- strata.tmp$coding
        ## old code:
        ## inla.eval(paste("strata.tmp = inla.strata(data.f$", strata.var, ")", sep=""))
        ## inla.eval(paste("data.f$", strata.var, " = strata.tmp$strata"))
        ## inla.eval(paste("data.l$baseline.hazard.strata.coding = strata.tmp$coding"))
    }

    ## if -1 the intercept is not included
    intercept <- inla.ifelse(attr(terms(formula), "intercept") == 0, FALSE, TRUE)

    ## we have to check if we have rw1/2,  or the iid
    if (inla.one.of(cont.hazard$model, c("rw1", "rw2"))) {
        d <- inla.ifelse((is.null(cont.hazard$diagonal) && cont.hazard$constr),
                         inla.set.f.default()$diagonal, cont.hazard$diagonal)
        sc <- paste0(", scale.model = ", inla.ifelse(
                                             is.null(cont.hazard$scale.model),
                                             inla.getOption("scale.model.default"),
                                             cont.hazard$scale.model
                                         ))
    } else {
        sc <- ""
        d <- if (is.null(cont.hazard$diagonal)) 0.0 else cont.hazard$diagonal
    }

    f.hazard <- paste(
        "~",
        inla.ifelse(intercept, "1 +", "-1 +"),
        ". + f(baseline.hazard, model=\"", cont.hazard$model, "\"",
        ", values = baseline.hazard.values",
        ", hyper = list(prec=", as.character(cont.hazard$hyper), ")",
        ", constr = ", cont.hazard$constr,
        ", diagonal = ", d,
        sc,
        inla.ifelse(is.null(strata.var), "", paste(", replicate=", strata.var)),
        ")",
        sep = ""
    )

    new.formula <- update(update(formula, as.formula(f.hazard)), y..coxph ~ .)

    data.list <- c(res$data.list, data.l)
    if (!is.null(strata.tmp)) {
        data.list <- c(data.list, list(baseline.hazard.strata.coding = strata.tmp$coding))
    }

    stopifnot(is.character(tag))
    if (nchar(tag) > 0) {
        old.names <- c(
            "y..coxph",
            "expand..coxph",
            "baseline.hazard", 
            "baseline.hazard.idx", 
            "baseline.hazard.time", 
            "baseline.hazard.length"
        )
        new.names <- c(
            paste0("y", tag, "..coxph"),
            paste0("expand", tag, "..coxph"),
            paste0("baseline", tag, ".hazard"), 
            paste0("baseline", tag, ".hazard.idx"), 
            paste0("baseline", tag, ".hazard.time"), 
            paste0("baseline", tag, ".hazard.length")
        )
        ## formula
        nf.text <- as.character(new.formula)
        nf.text <- paste(nf.text[2], nf.text[1], nf.text[3])
        for (i in seq_along(old.names)) {
            nf.text <- gsub(old.names[i], new.names[i], nf.text, fixed = TRUE)
        }
        new.formula <- as.formula(nf.text)

        ## data
        rd.names <- names(res$data)
        for (i in seq_along(old.names)) {
            rd.names <- gsub(old.names[i], new.names[i], rd.names, fixed = TRUE)
        }
        names(res$data) <- rd.names

        ## data.list
        dl.names <- names(data.list)
        for (i in seq_along(old.names)) {
            dl.names <- gsub(old.names[i], new.names[i], dl.names, fixed = TRUE)
        }
        names(data.list) <- dl.names
    }

    expand..coxph.idx <- which(names(res$data) == paste0("expand", tag, "..coxph"))

    return(list(
        formula = new.formula,
        data = res$data,
        data.list = data.list,
        family = "poisson",
        E = res$data[, "E..coxph"],
        expand.df = res$data[, expand..coxph.idx],
        control.hazard = cont.hazard
    ))
}


#' @rdname coxph
#' @export inla.rbind.data.frames
`inla.rbind.data.frames` <- function(...)
{
    ## rbind data.frames padding with NA

    if (inla.require("dplyr")) {
        ## a better implementation is \code{dplyr::bind_rows})
        return(dplyr::bind_rows(...))
    } else {
        rbind.two.data.frames <- function(df1, df2) {
            df <- NULL
            nr1 <- nrow(df1)
            nr2 <- nrow(df2)
            nams <- sort(unique(c(names(df1), names(df2))), decreasing = TRUE)

            for (nam in nams) {
                if (nam %in% names(df1)) {
                    idx <- which(names(df1) %in% nam)
                    if (nam %in% names(df)) {
                        idxx <- which(names(df) %in% nam)
                        df[1:nr1, idxx] <- df1[1:nr1, idx]
                    } else {
                        if (!is.null(df)) {
                            df <- data.frame(c(df1[1:nr1, idx], rep(NA, nr2)), df)
                        } else {
                            df <- data.frame(c(df1[1:nr1, idx], rep(NA, nr2)))
                        }
                        names(df)[1] <- nam
                    }
                }
                if (nam %in% names(df2)) {
                    idx <- which(names(df2) %in% nam)
                    if (nam %in% names(df)) {
                        idxx <- which(names(df) %in% nam)
                        df[nr1 + 1:nr2, idxx] <- df2[1:nr2, idx]
                    } else {
                        if (!is.null(df)) {
                            df <- data.frame(c(rep(NA, nr1), df2[1:nr2, idx]), df)
                        } else {
                            df <- data.frame(c(rep(NA, nr1), df2[1:nr2, idx]))
                        }
                        names(df)[1] <- nam
                    }
                }
            }

            return(df)
        }

        args <- list(...)
        stopifnot(length(args) >= 2L)

        for (i in 2L:length(args)) {
            if (i == 2L) {
                df <- rbind.two.data.frames(args[[1]], args[[2]])
            } else {
                df <- rbind.two.data.frames(df, args[[i]])
            }
        }

        return(df)
    }
}
