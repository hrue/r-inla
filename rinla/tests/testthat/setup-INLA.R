testthat_eval_env <- environment()
testthat_eval_env <- tryCatch(
    teardown_env(),
    error = function(e) { testthat_eval_env }
)

withr::local_options(
    lifecycle_verbosity = "warning",
    .local_envir = testthat_eval_env
)

local_inla_options(
    num.threads = "1:1",
    fmesher.evolution.warn = TRUE,
    fmesher.evolution.verbosity = "warn",
    .envir = testthat_eval_env
)
