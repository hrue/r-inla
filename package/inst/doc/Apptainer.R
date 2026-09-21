## ----apptainer check, eval=FALSE----------------------------------------------
# $ apptainer version

## ----apptainer, eval=FALSE----------------------------------------------------
# # To pull the docker image
# $ apptainer pull your_container.sif docker://rocker/geospatial:latest
# # once the sif is downloaded, get into an interactive shell
# $ apptainer shell your_container.sif
# # Inside the interactive shell, one can install INLA on a personal library path
# $ R --verbose
# # Now, in a R environment
# > options(repos = c(
#     INLA = 'https://inla.r-inla-download.org/R/testing',
#     CRAN = 'https://cloud.r-project.org'))
# > install.packages("INLA")
# # One may be asked to create a personal library path.
# > install.packages("inlabru")
# # quit R
# > q()

## ----exec, eval=FALSE---------------------------------------------------------
# # To execute an R file with apptainer
# $ apptainer exec ./your_container.sif Rscript --no-restore --no-save --verbose file_to_run.R

