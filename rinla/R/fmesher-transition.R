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
    
    w <- isTRUE(inla.getOption("fmesher.evolution.warn"))
    verb <- inla.getOption("fmesher.evolution.verbosity")
    
    verb <- match.arg(verb, c("default", "soft", "warn", "stop"))
    level <- match.arg(level, c("default", "soft", "warn", "stop"))

    if (identical(level, "stop") || (w && identical(verb, "stop"))) {
        lifecycle::deprecate_stop(
            when = when,
            what = what,
            with = with,
            details = details,
            env = env
        )
    } else if (w && (identical(level, "warn") || identical(verb, "soft"))) {
        lifecycle::deprecate_warn(
            when = when,
            what = what,
            with = with,
            details = details,
            id = id,
            always = always,
            env = env,
            user_env = user_env
        )
    } else if (w && (identical(level, "soft") || identical(verb, "soft"))) {
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
