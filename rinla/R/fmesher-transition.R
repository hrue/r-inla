#' @title Helper for locally setting inla options
#' 
#' @description Locally sets INLA options with [inla.setOption()] and restores
#' when the calling function terminates, like [withr::local_options()] does for
#' global R options.
#' @param \dots Named arguments to be passed to [INLA::inla.setOption()].
#' @param .envir The environment in which the options should be set. Defaults to
#'   the calling environment.
#' @param .save_only If `TRUE`, the options are not set, but the current values
#'   are saved and will be restored when the calling function terminates, even if
#'   other code sets them directly by calling [inla.setOption()]. This can be
#'   used to protect specific options against local changes that should not have
#'   global effect. Default is `FALSE`, meaning that the options are set.
#' 
#' @returns logical; If the option setting was successful, `TRUE` is returned,
#' otherwise `FALSE`.
#' @examples
#' print(paste0("Before: ", INLA::inla.getOption("num.threads")))
#' local({
#'   fun1 <- function() { INLA::inla.setOption("num.threads" = "1:1") }
#'   fun2 <- function() {
#'     local_inla_options(num.threads = NULL, .save_only = TRUE)
#'     print(paste0("Local value: ", INLA::inla.getOption("num.threads")))
#'     fun1()
#'     print(paste0("Local changed value: ",
#'                  INLA::inla.getOption("num.threads")))
#'     invisible()
#'   }
#'   
#'   fun2()
#'   print(paste0("Restored: ", INLA::inla.getOption("num.threads")))
#'   local_inla_options(num.threads = "2:1")
#'   print(paste0("Change: ", INLA::inla.getOption("num.threads")))
#'   fun2()
#'   print(paste0("Restored: ",
#'                INLA::inla.getOption("num.threads")))
#' })
#' print(paste0("After: ", INLA::inla.getOption("num.threads")))
#' @export
local_inla_options <- function(...,
                               .envir = parent.frame(),
                               .save_only = FALSE) {
    # Set INLA options
    inla_options <- list(...)
    old_inla_options <- list()
    if (length(inla_options) > 0) {
        for (name in names(inla_options)) {
            old_inla_options[[name]] <- tryCatch(
                INLA::inla.getOption(name),
                error = function(e) {
                    e
                }
            )
            
            if (inherits(old_inla_options[[name]], "simpleError")) {
                return(FALSE)
            }
            
            if (!.save_only) {
                e <- tryCatch(
                    INLA::inla.setOption(name, inla_options[[name]]),
                    error = function(e) {
                        e
                    }
                )
                if (inherits(e, "simpleError")) {
                    return(FALSE)
                }
            }
            
            withr::defer(
                INLA::inla.setOption(name, old_inla_options[[name]]),
                .envir
            )
        }
    }
    TRUE
}


## Helpers to avoid a runtime dependency on testthat.
## See documentation ?testthat::is_testing
testthat_is_testing <- function() {
    identical(Sys.getenv("TESTTHAT"), "true")
}
testthat_is_snapshot <- function() {
    identical(Sys.getenv("TESTTHAT_IS_SNAPSHOT"), "true")
}
testthat_is_checking <- function() 
{
    identical(Sys.getenv("TESTTHAT_IS_CHECKING"), "true")
}
testthat_is_parallel <- function() 
{
    identical(Sys.getenv("TESTTHAT_IS_PARALLEL"), "true")
}
testthat_testing_package <- function() 
{
    Sys.getenv("TESTTHAT_PKG")
}

fmesher_deprecate <- function(level = NULL,
                              evo = NULL,
                              when,
                              what,
                              with = NULL,
                              details = NULL,
                              id = NULL,
                              always = FALSE,
                              env = rlang::caller_env(),
                              user_env = rlang::caller_env(2)) {
    w <- isTRUE(inla.getOption("fmesher.evolution.warn")) ||
        (testthat_is_testing() && !testthat_is_snapshot())
    verb <- inla.getOption("fmesher.evolution.verbosity")
    
    verb <- match.arg(verb, c("default", "soft", "warn", "stop"))
    level <- match.arg(level, c("default", "soft", "warn", "stop"))
    
    details <- c(
        details,
        "For more information, see https://inlabru-org.github.io/fmesher/articles/inla_conversion.html",
        "To silence these deprecation messages in old legacy code, set `inla.setOption(fmesher.evolution.warn = FALSE)`.",
        "To ensure visibility of these messages in package tests, also set `inla.setOption(fmesher.evolution.verbosity = 'warn')`."
    )
    
    if (w) {
        if (identical(verb, "default")) {
            verb <- "soft"
        }
        if (identical(verb, "soft") && testthat_is_testing()) {
            verb <- "warn"
        }
        withr::local_options(
            lifecycle_verbosity = c(
                soft = "default",
                warn = "warning",
                stop = "error"
            )[verb]
        )
    } else {
        withr::local_options(lifecycle_verbosity = "quiet")
        if (identical(verb, "default")) {
            verb <- "soft"
        }
    }
    if (identical(level, "default")) {
        level <- verb
    }
    
    if (identical(level, "stop")) {
        lifecycle::deprecate_stop(
            when = when,
            what = what,
            with = with,
            details = details,
            env = env
        )
    } else if (identical(level, "warn")) {
        lifecycle::deprecate_warn(
            when = when,
            what = what,
            with = with,
            details = details,
            id = id,
            always = always ||
                (testthat_is_testing() &&
                     !testthat_is_snapshot()),
            env = env,
            user_env = user_env
        )
    } else if (identical(level, "soft")) {
        lifecycle::deprecate_soft(
            when = when,
            what = what,
            with = with,
            details = details,
            id = id,
            env = env,
            user_env = user_env
        )
    }

    return(TRUE)
}
