`inla.version` = function(details=FALSE, quiet =FALSE, hg=FALSE) {
    if (hg && missing(details) && missing(quiet))
        return ("hgid not available.")
    cat("\n",
        "This version of R-INLA use the code from inla.googlecode.com",
        "\n",
        "No version-number is available.",
        "\n",
        "\n")
    return (invisible())
}
