fmesher_deprecate_soft <- function(evo = NULL,
                                   when,
                                   what,
                                   with = NULL,
                                   details = NULL,
                                   id = NULL,
                                   env = rlang::caller_env(),
                                   user_env = rlang::caller_env(2)) {
    if (!is.null(evo) && (inla.getOption("fmesher.evolution") < evo)) {
        return(FALSE)
    }
    if (inla.getOption("fmesher.evolution.warn")) {
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
fmesher_deprecate_warn <- function(evo = NULL,
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
    if (inla.getOption("fmesher.evolution.warn")) {
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
    }
    return(TRUE)
}
fmesher_deprecate_stop <- function(evo = NULL,
                                   when,
                                   what,
                                   with = NULL,
                                   details = NULL,
                                   env = rlang::caller_env()) {
    if (!is.null(evo) && (inla.getOption("fmesher.evolution") < evo)) {
        return(FALSE)
    }
    lifecycle::deprecate_stop(
        when = when,
        what = what,
        with = with,
        details = details,
        env = env
    )
}
