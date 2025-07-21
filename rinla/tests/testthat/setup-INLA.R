# Save the num.threads option so it can be restored
testthat_inla_num_threads <- INLA::inla.getOption("num.threads")
# Set num.threads = "1:1" to ensure within-system repeatability
INLA::inla.setOption(num.threads = "1:1")

testthat_fmesher_evolution_warn <- INLA::inla.getOption("fmesher.evolution.warn")
testthat_fmesher_evolution_verb <- INLA::inla.getOption("fmesher.evolution.verbosity")

INLA::inla.setOption(fmesher.evolution.warn = TRUE)
INLA::inla.setOption(fmesher.evolution.verbosity = "warn")

testthat_eval_env <- environment()
testthat_eval_env <- tryCatch(teardown_env(), error = function(e) testthat_eval_env)

withr::local_options(lifecycle_verbosity = "warning", .local_envir = testthat_eval_env)
withr::defer(
    INLA::inla.setOption(
        num.threads = testthat_inla_num_threads,
        fmesher.evolution.warn = testthat_fmesher_evolution_warn,
        fmesher.evolution.verbosity = testthat_fmesher_evolution_verb
    ),
    testthat_eval_env
)
