# Save the num.threads option so it can be restored
testthat_inla_num_threads <- INLA::inla.getOption("num.threads")
# Set num.threads = "1:1" to ensure within-system repeatability
INLA::inla.setOption(num.threads = "1:1")

withr::defer(
    INLA::inla.setOption(num.threads = testthat_inla_num_threads),
    teardown_env()
)

testthat_fmesher_evolution_warn <- INLA::inla.getOption("fmesher.evolution.warn")
testthat_fmesher_evolution_verb <- INLA::inla.getOption("fmesher.evolution.verbosity")

INLA::inla.setOption(fmesher.evolution.warn = TRUE)
withr::defer(
    INLA::inla.setOption(fmesher.evolution.warn = testthat_fmesher_evolution_warn),
    teardown_env()
)
INLA::inla.setOption(fmesher.evolution.verbosity = "soft")
withr::defer(
    INLA::inla.setOption(fmesher.evolution.verbosity = testthat_fmesher_evolution_verb),
    teardown_env()
)
