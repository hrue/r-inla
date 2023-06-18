# Save the num.threads option so it can be restored
testthat_inla_num_threads <- INLA::inla.getOption("num.threads")
# Set num.threads = "1:1" to ensure within-system repeatability
INLA::inla.setOption(num.threads = "1:1")

withr::defer(
    INLA::inla.setOption(num.threads = testthat_inla_num_threads),
    teardown_env()
)
