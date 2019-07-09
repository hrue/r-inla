## Export: inla.pardiso inla.pardiso.check

##!\name{inla.pardiso}
##!\alias{inla.pardiso}
##!\alias{inla.pardiso.check}
##!\title{PARDISO support in R-INLA}
##!\description{
##!Describe and check the PARDISO support in R-INLA
##!}
##!\usage{
##!inla.pardiso()
##!inla.pardiso.check()
##!}
##!\arguments{
##!}
##!\details{
##!\code{inla.pardiso()} describes the \code{PARDISO} support in R-INLA, how to get the license
##!key and enable it in the \code{R-INLA} package. \code{inla.pardiso.check()} check if the
##!\code{PARDISO} support is working.
##!}
##!\value{
##!}
##!\author{Havard Rue \email{hrue@r-inla.org}}

`inla.pardiso` = function()
{
    browseURL("https://pardiso-project.org/r-inla/")
    cat("\n\n",
        "The R-INLA page at the PARDISO project should open in your browser, \n",
        "which gives more information (https://pardiso-project.org/r-inla)\n",
        "\n", 
        "The PARDISO package is a thread-safe, high-performance, and easy\n",
        "to use software for solving sparse symmetric that arises in the\n",
        "R-INLA approach to Bayesian inference within the R project for \n",
        "Statistical Computing.\n",
        "\n",
        "R-INLA supports PARDISO for MacOSX and Linux.\n",
        "Windows support will be added later.\n", 
        "\n", 
        "To improve R-INLA within parallel and high-performance computing, \n",
        "we need a modern sparse-matrix library to build upon, and\n",
        "and for this, PARDISO, is the very best choice in our opinion.\n",
        "\n",
        "The PARDISO library is not enabled nor loaded by default\n",
        "when calling inla(). You can you can enable it, by\n",
        "\n", 
        "1. Obtaining a license key from pardiso-project.org/r-inla\n",
        "   The licenses is tied up to the username, so the same license\n",
        "   can be used one all computers where you have the same username.\n",
        "   (You can also collect various licenses in the same file, each\n",
        "   connected with different usernames.)\n", 
        "2. Add the license key to R-INLA, by assigning the full path\n",
        "   of the key-file to variable 'pardiso.license', for example,\n", 
        "     library(INLA)\n",
        "     inla.setOption(pardiso.license=\"~/sys/licenses/pardiso.lic\")\n",
        "3. That is it!\n", 
        "4. You can check PARDISO setup, by\n",
        "     inla.pardiso.check()\n", 
        "\n",
        "Using remote computing with inla.call=\"remote\", you may need to add\n",
        "the full path to the license file (on the remote host) in the\n",
        "~/.inlarc file:\n",
        "     PardisoLicenseFile=\"FullPathToTheLicenseFile\"\n",
        "Otherwise, your local license file will be used if its available.", 
        "\n", 
        "By default, the PARDISO solver will run in 'serial' mode, but for\n",
        "large models, the parallel mode is preferred.\n",
        "The mode can be controlled with\n", 
        "     control.compute=list(openmp.strategy=\"pardiso.serial\"), \n", 
        "and\n", 
        "     control.compute=list(openmp.strategy=\"pardiso.parallel\"), \n",
        "Nested parallelism, will be supported in the future.\n", 
        "\n",
        "WARNING: The PARDISO support is highly experimental for the moment\n",
        "and changes in the interface will likely occour.\n",
        "\n",
        "Please report problems to help@r-inla.org\n",
        "\n", 
        "Havard Rue & Olaf Schenk\n",
        "May 2018\n",
        "\n")
}

`inla.pardiso.check` = function() 
{
    t.dir = inla.tempdir()
    inla.set.sparselib.env(inla.dir = t.dir)
    if (inla.os("linux") || inla.os("mac")) {
        ret = system(paste(shQuote(inla.getOption("inla.call")), "-m pardiso"), intern=TRUE)
    } else if(inla.os("windows")) {
        ret = system(paste(shQuote(inla.getOption("inla.call")), "-m pardiso"), intern=TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }
    unlink(t.dir, recursive = TRUE)

    return (ret)
}
