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

fmesher_deprecate_allow <- function(evo = NULL,
                                    env = rlang::caller_env(),
                                    user_env = rlang::caller_env(2)) {
    if (!is.null(evo) && (inla.getOption("fmesher.evolution") < evo)) {
        return(FALSE)
    }
    return(TRUE)
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
    if (!is.null(evo) && (inla.getOption("fmesher.evolution") < evo)) {
        return(FALSE)
    }
    
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
