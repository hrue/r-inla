source("../../../rinla/R/models.R")
source("../../../rinla/R/utils.R")
m = inla.models()
for(nm in names(m)) {
    if (!file.exists(nm))
        dir.create(nm)
    for (mod  in names(m[[nm]])) {
        fnm = paste(nm, "/", mod, ".tex",  sep="")
        cat("Generate file", fnm, "\n")
        inla.models.generate.tex(m[[nm]][[mod]], fnm)
    }
}
