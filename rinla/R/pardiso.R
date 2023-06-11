#' @title Describe and check the PARDISO support in R-INLA
#' 
#' @description
#' `inla.pardiso()` describes the `PARDISO` support in R-INLA, how to
#' get the license key and enable it in the `R-INLA` package.
#' `inla.pardiso.check()` check if the `PARDISO` support is working.
#' 
#' @name pardiso
#' @aliases pardiso inla.pardiso inla.pardiso.check
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname pardiso
#' @export inla.pardiso

`inla.pardiso` <- function() {
    stop("The use of the PARDISO library in R-INLA is currently paused.\n       Old license keys can be used until they expire.")
}

`inla.pardiso.INVALID.DESCRIPTION` <- function() {
    browseURL("https://pardiso-project.org/r-inla/")
    cat(
        "\n\n",
        "The R-INLA page at the PARDISO project should open in your browser, \n",
        "which gives more information (https://pardiso-project.org/r-inla)\n",
        "\n",
        "The PARDISO package is a thread-safe, high-performance, and easy\n",
        "to use software for solving sparse symmetric that arises in the\n",
        "R-INLA approach to Bayesian inference within the R project for \n",
        "Statistical Computing.\n",
        "R-INLA supports PARDISO for MacOSX and Linux.\n",
        "Windows is not supported for the moment.\n",
        "\n",
        "To improve R-INLA within parallel and high-performance computing, \n",
        "we need a modern sparse-matrix library to build upon, and\n",
        "and for this, PARDISO, is the very best choice in our opinion.\n",
        "\n",
        "The PARDISO library is not enabled by default\n",
        "when calling inla(). You can you can enable it, by\n",
        "\n",
        "1. Obtaining a license key from pardiso-project.org/r-inla\n",
        "   The licenses is tied up to the username, so the same license\n",
        "   can be used one all computers where you have the same username.\n",
        "   (You can also collect various licenses in the same file, each\n",
        "   connected with different usernames.)\n",
        "2. Add the license key to R-INLA, by assigning the full path\n",
        "   of the license-file to variable 'pardiso.license', for example,\n",
        "     library(INLA)\n",
        "     inla.setOption(pardiso.license=\"~/sys/licenses/pardiso.lic\")\n",
        "   (Alternatively, the key(s) can be added directly,  like\n",
        "     inla.setOption(pardiso.license=\"303BCDF0C08EF46A51...A81DA\")\n",
        "   Use newline to separate more than one key.)\n", 
        "3. That is it!\n",
        "4. You can check that PARDISO is enabled and is working, by\n",
        "     inla.pardiso.check()\n",
        "\n",
        "Using remote computing with inla.call=\"remote\", you may need to add\n",
        "the full path to the license file (on the remote host) in the\n",
        "~/.inlarc file:\n",
        "     PardisoLicenseFile=\"FullPathToTheLicenseFile\"\n",
        "Otherwise, your local license file will be used if its available.\n",
        "\n",
        "You can also request PARDISO explicitely with \n",
        "     control.compute=list(smtp=\"pardiso\"), \n",
        "Nested parallelism, is controlled by the \"num.threads\" argument\n",
        "     num.threads=\"A:B\" \n",
        "with A threads in the outer layer and B threads in the inner.\n",
        "(Think of this as dealing with A matrices at the same time, and \n",
        "using B threads on each of them.)\n",
        "Good choices are\n",
        "-  A=number.of.hyperparameters+1 or A=2*number.of.hyperparameters\n",
        "-  B=a low number like 1, 2, 4, ... higher for huge models.\n",
        "You also need to take the number of cores/threads available into account.\n",
        "A zero can be used for both A and B and their values will be set by the program.\n",
        "\n",
        "PARDISO support is still under development, and \n",
        "changes in the interface/behaviour might occour.\n",
        "\n",
        "Please report problems to help@r-inla.org\n",
        "\n",
        "Havard Rue & Olaf Schenk\n",
        "May 2018\n",
        "Updated August 2020\n",
        "Updated March 2022\n",
        "\n"
    )
}

#' @rdname pardiso
#' @export inla.pardiso.check
`inla.pardiso.check` <- function() {
    t.dir <- inla.tempdir()
    inla.set.sparselib.env(inla.dir = t.dir)
    if (inla.os("linux") || inla.os("mac") || inla.os("mac.arm64")) {
        ret <- system(paste(shQuote(inla.call.no.remote()), "-m pardiso"), intern = TRUE)
    } else if (inla.os("windows")) {
        ret <- system(paste(shQuote(inla.call.no.remote()), "-m pardiso"), intern = TRUE)
    } else {
        stop("\n\tNot supported architecture.")
    }
    unlink(t.dir, recursive = TRUE)

    return(ret)
}
