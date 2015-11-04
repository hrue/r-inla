options(keep.source=TRUE)
source("../../../rinla/R/inlaEnv.R")
source("../../../rinla/R/models.R")
source("../../../rinla/R/models-generate.R")
source("../../../rinla/R/utils.R")
source("../../../rinla/R/version.R")
m = inla.models()
for(nm in names(m)) {
    if (!file.exists(nm))
        dir.create(nm)
    for (mod  in names(m[[nm]])) {
        fnm = paste(nm, "/", mod, ".tex",  sep="")
        ## remove ':' in the filename
        fnm = gsub(":", "", fnm)
        cat("Generate file", fnm, "\n")
        inla.models.generate.tex(m[[nm]][[mod]], fnm)
    }
}
