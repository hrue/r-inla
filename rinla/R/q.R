#' @title Control and view a remote inla-queue
#' 
#' @description
#' Control and view a remote inla-queue of submitted jobs
#' 
#' `inla.qstat` show job(s) on the server, `inla.qget` fetch the
#' results (and by default remove the files on the server), `inla.qdel`
#' removes a job on the server and `inla.qnuke` remove all jobs on the
#' server. `inla.qlog` fetches the logfile only.
#' 
#' The recommended procedure is to use `r=inla(..., inla.call="submit")`
#' and then do `r=inla.qget(r)` at a later stage.  If the job is not
#' finished, then `r` will not be overwritten and this step can be
#' repeated.  The reason for this procedure, is that some information usually
#' stored in the result object does not go through the remote server, hence
#' have to be appended to the results that are retrieved from the server. Hence
#' doing `r=inla(..., inla.call="submit")` and then later retrive it using
#' `r=inla.qget(1)`, say, then `r` does not contain all the usual
#' information.  All the main results are there, but administrative information
#' which is required to call `inla.hyperpar` or `inla.rerun` are not
#' there.
#' 
#' @name inla.qstat
#' @aliases inla.qstat inla.q inla.qget inla.qdel inla.qnuke inla.qlog
#' summary.inla.q print.inla.q
#' @param id The job-id which is the output from `inla` when the job is
#' submitted, the job-number or job-name. For `inla.qstat`, `id` is
#' optional and if omitted all the jobs will be listed.
#' @param remove Logical If FALSE, leave the job on the server after retrival,
#' otherwise remove it (default).
#' @param x An `inla.q`-object which is the output from `inla.qstat`
#' @param object An `inla.q`-object which is the output from
#' `inla.qstat`
#' @param ...  other arguments.
#' @return `inla.qstat` returns an `inla.q`-object with information
#' about current jobs.
#' @author Havard Rue
#' @seealso [inla()]
#' @examples
#' 
#' \dontrun{
#' r = inla(y~1, data = data.frame(y=rnorm(10)), inla.call="submit")
#' inla.qstat()
#' r = inla.qget(r, remove=FALSE)
#' inla.qdel(1)
#' inla.qnuke()
#' }
#' 
#' @rdname q
NULL

#' @rdname q
#' @method summary inla.q
#' @export
`summary.inla.q` <- function(object, ...) {
    print(object, ...)
}

#' @rdname q
#' @method print inla.q
#' @export
`print.inla.q` <- function(x, ...) {
    for (k in seq_along(x)) {
        cat("Job:", x[[k]]$no, "  Id:", x[[k]]$id, "  Size:", x[[k]]$size, "  Status:", x[[k]]$status, "\n")
    }
    return(invisible(x))
}

#' @rdname q
#' @export
`inla.qget` <- function(id, remove = TRUE) {
    return(inla.q(cmd = "get", id = id, remove = remove))
}

#' @rdname q
#' @export
`inla.qdel` <- function(id) {
    return(inla.q(cmd = "del", id = id))
}



#' @rdname q
#' @export
`inla.qstat` <- function(id) {
    return(inla.q(cmd = "stat", id = id))
}

#' @rdname q
#' @export
`inla.qlog` <- function(id) {
    return(inla.q(cmd = "log", id = id))
}

#' @rdname q
#' @export
`inla.qnuke` <- function() {
    return(inla.q(cmd = "nuke"))
}

`inla.q` <- function(cmd = c("get", "del", "stat", "log", "nuke"), id, remove = TRUE) {
    cmd <- match.arg(cmd)

    ## do a quick return here if possible
    if (cmd == "log" && missing(id)) {
          return(NULL)
      }

    ## define some environment variables for remote computing
    vars <- c(
        INLA_PATH = system.file("bin", package = "INLA"),
        INLA_OS = inla.os.type(),
        INLA_VERSION = inla.version("version")
    )
    if (inla.os("windows")) {
        vars <- c(vars,
            INLA_SSH_AUTH_SOCK = inla.getOption("ssh.auth.sock"),
            INLA_CYGWIN_HOME = inla.getOption("cygwin.home"),
            INLA_HOME = inla.cygwin.map.filename(gsub("\\\\", "/", inla.get.HOME()))
        )
    } else {
        vars <- c(vars,
            INLA_HOME = inla.get.HOME()
        )
        if (Sys.getenv("SSH_AUTH_SOCK") == "") {
            vars <- c(vars,
                INLA_SSH_AUTH_SOCK = inla.getOption("ssh.auth.sock")
            )
        }
    }
    do.call("Sys.setenv", as.list(vars))

    inla.call <- system.file("bin/remote/inla.q", package = "INLA")
    if (inla.os("windows")) {
        inla.call <- paste(inla.call, ".cygwin", sep = "")
    }

    if (!missing(id) && is.list(id)) {
        ret <- id$ret
        id <- id$id
    } else {
        ret <- NULL
    }

    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        output <- system(paste(shQuote(inla.call), cmd, inla.ifelse(missing(id) || is.null(id), "NULL", id), as.numeric(remove)), intern = TRUE)
    } else {
        bat.file <- paste(tempfile(), ".BAT", sep = "")
        cat("@echo off\n", file = bat.file, append = FALSE)
        cat(paste(shQuote(inla.call), cmd, inla.ifelse(missing(id) || is.null(id), "NULL", id), as.numeric(remove)))
        output <- try(shell(paste("@", shQuote(bat.file), sep = ""), wait = TRUE, intern = TRUE))
        unlink(bat.file)
    }

    if (length(grep("^ERROR", output)) > 0) {
        output <- gsub("^ERROR *", "", output)
        stop(output)
    } else if (length(grep("^DELETE", output)) > 0) {
        output <- gsub("^DELETE *", "", output)
        cat("Delete job", output, "\n")
    } else if (length(grep("^NUKE", output)) > 0) {
        output <- gsub("^NUKE *", "", output)
    } else if (length(grep("^LOG", output)) > 0) {
        output <- gsub("^LOG *", "", output)
        if (file.exists(output)) {
            logfile <- list(logfile = gsub("\t", "        ", readLines(output)))
            try(unlink(output, recursive = TRUE), silent = TRUE)
        } else {
            logfile <- NULL
        }
        return(logfile)
    } else if (length(output) > 0 && length(strsplit(output, " ")[[1]]) == 1) {
        r <- inla.collect.results(output, file.log = paste(output, "/results.files/Logfile.txt", sep = ""))
        rr <- c(r, ret)
        class(rr) <- class(r)
        if (!is.null(ret) && ret$.args$keep == TRUE && !is.null(ret$.args$working.directory)) {
            ## copy the files to the target directory
            d.fnm <- paste(ret$.args$working.directory, "/results.files", sep = "")
            inla.dir.create(d.fnm)
            files.to.copy <- paste(output, "/", dir(output, recursive = TRUE), sep = "")
            file.copy(files.to.copy, d.fnm, recursive = TRUE)
        } else {
            try(unlink(output, recursive = TRUE), silent = TRUE)
        }
        if (!is.null(ret) && ret$.args$keep == FALSE && is.null(ret$.args$working.directory)) {
            try(unlink(ret$misc$inla.dir, recursive = TRUE), silent = TRUE)
        }
        return(rr)
    } else {
        if (length(output) >= 1 && nchar(output[1]) > 0) {
            output <- lapply(
                strsplit(output, " +"),
                function(a) {
                    names(a) <- c("id", "no", "pid", "status", "size")
                    return(as.list(a))
                }
            )
        } else {
            output <- list()
        }
        class(output) <- "inla.q"
        return(output)
    }
}
