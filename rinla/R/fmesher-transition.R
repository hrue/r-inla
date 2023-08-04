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
    w <- inla.getOption("fmesher.evolution.warn")
    if (is.logical(w)) {
        if (w) {
            w <- "soft"
        } else {
            w <- "none"
        }
    }
    w <- match.arg(w, c("none", "soft", "warn", "stop"))
    level <- match.arg(level, c("none", "soft", "warn", "stop"))
    levels <- c(level, w)

    if (!is.null(evo) && (inla.getOption("fmesher.evolution") < evo)) {
        return(FALSE)
    }

    if ("stop" %in% levels) {
        lifecycle::deprecate_stop(
            when = when,
            what = what,
            with = with,
            details = details,
            env = env
        )
    } else if ("warn" %in% levels) {
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
    } else if ("soft" %in% levels) {
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
