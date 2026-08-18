#' @title Add sTiles support in R-INLA
#' 
#' @description
#' `inla.stiles()` opens the web-page with
#' description how to enable `stiles` support in R-INLA
#' 
#' @name stiles
#' @aliases stiles inla.stiles
#' @author Havard Rue \email{hrue@@r-inla.org}
#' @rdname stiles
#' @export inla.stiles

`inla.stiles` <- function() {
    browseURL("https://github.com/hrue/r-inla/releases/latest/download/")
    cat('\n',"
1. Download the matching binary file and store it in
   directory/folder <D>. (You have to replace <D> with
   the actual path.)
   - https://github.com/hrue/r-inla/releases/latest/download/inla-linux-x86_64-portable.tar.gz
   - https://github.com/hrue/r-inla/releases/latest/download/inla-linux-x86_64-v3-portable.tar.gz
   - https://github.com/hrue/r-inla/releases/latest/download/inla-linux-arm64-armv82-portable.tar.gz
   - https://github.com/hrue/r-inla/releases/latest/download/inla-macos-arm64-portable.tar.gz
   - https://github.com/hrue/r-inla/releases/latest/download/inla-macos-x86_64-portable.tar.gz
   - https://github.com/hrue/r-inla/releases/latest/download/inla-windows-x86_64.zip
2. Go to <D> and un-pack/zip the file
3. Add this to the top of your R-script:
       inla.setOption(inla.call='<D>/bin/inla', smtp='stiles')     ## Mac & Linux
       inla.setOption(inla.call='<D>/bin/inla.exe', smtp='stiles') ## Windows
4. Next time you upgrade R-INLA you have to do similar.\n\n")

    ans <- ""
    sysname <- Sys.info()["sysname"]
    machine <- Sys.info()["machine"]
    if (sysname == "Linux") {
        if (machine == "x86_64") {
            s <- system("getconf GNU_LIBC_VERSION", intern = TRUE)
            s <- as.numeric(strsplit(s, " ")[[1]][2])
            newer <- (s >= 2.38)
            if (newer) {
                ans <- "inla-linux-x86_64-v3-portable"
            } else {
                ans <- "inla-linux-x86_64-portable"
            }
        } else if (machine== "aarch64" || machine == "arm64") {
            ans <- "inla-linux-arm64-armv82-portable"
        }
    } else if (sysname == "Darwin") {
        if (machine == "x86_64") {
            ans <- "inla-macos-x86_64-portable"
        } else {
            ans <- "inla-macos-arm64-portable"
        }
    } else if (sysname == "Windows") {
        ans <- "inla-windows-x86_64"
    }
    if (ans != "") {
        cat("\n==> I think '", ans, "' is your choice\n", sep = "")
    }

    cat("\n==> More information is available here: ", 
        "https://github.com/hrue/r-inla/releases/tag/v", inla.version("version"), "\n", sep = "")

    return (invisible())
}
