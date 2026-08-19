#' @title Install a pre-built INLA binary with sTiles support
#'
#' @description
#' `inla.stiles.install()` performs what [inla.stiles()] describes: it works
#' out which pre-built binary matches this machine, downloads it from the
#' releases page, unpacks it into a cache directory, checks that it runs, and
#' points `inla.call` at it. One call instead of the manual download,
#' un-pack and `inla.setOption()` steps.
#'
#' The binaries are portable: six of them cover Linux (x86-64 baseline and
#' x86-64-v3), Linux arm64, macOS (Intel and Apple silicon) and Windows, so the
#' choice follows from the operating system, the architecture, the C library
#' version and -- on x86-64 -- the instruction sets the processor actually
#' supports. Each bundle embeds a matching `libstiles`, and carries a
#' `BUILDINFO` file recording the compiler, flags and library versions that
#' produced it.
#'
#' @param tag     Release to install from, e.g. `"v26.08.18"`. Defaults to the
#'                current release. Pin it for reproducibility.
#' @param dir     Where to unpack. Defaults to a per-user cache directory.
#' @param force   Re-download even when the binary is already installed.
#' @param smtp    Also select the sTiles sparse-matrix backend (default `TRUE`).
#' @param repo    GitHub repository holding the releases.
#' @param verbose Report each step.
#'
#' @return The path of the installed binary, invisibly.
#'
#' @examples
#' \dontrun{
#' inla.stiles.install()                  # newest release
#' inla.stiles.install(tag = "v26.08.18") # a specific one
#' }
#'
#' @seealso [inla.stiles()], [inla.binary.install()]
#' @name stiles.install
#' @aliases inla.stiles.install
#' @rdname stiles.install
#' @export inla.stiles.install

`inla.stiles.install` <- function(tag = NULL,
                                  dir = NULL,
                                  force = FALSE,
                                  smtp = TRUE,
                                  repo = "hrue/r-inla",
                                  verbose = TRUE) {
    say <- function(...) if (verbose) cat("*", paste0(..., collapse = ""), "\n")

    sysname <- Sys.info()["sysname"]
    machine <- Sys.info()["machine"]
    arm <- machine %in% c("aarch64", "arm64")

    ## The C library sets which Linux bundle can run at all, and on x86-64 the
    ## processor sets which one is safe: x86-64-v3 needs AVX2/FMA/BMI, i.e.
    ## Haswell (2013) or newer. A recent distribution on an older processor is
    ## an ordinary combination, and picking the v3 build there would abort with
    ## an illegal instruction on the first vectorised kernel.
    glibc <- NA_real_
    v3 <- FALSE
    if (sysname == "Linux") {
        s <- tryCatch(system("getconf GNU_LIBC_VERSION", intern = TRUE),
                      error = function(e) NA_character_,
                      warning = function(w) NA_character_)
        if (!is.na(s[1])) glibc <- suppressWarnings(as.numeric(strsplit(s, " ")[[1]][2]))
        if (!arm && file.exists("/proc/cpuinfo")) {
            flags <- grep("^flags", readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)[1]
            v3 <- length(flags) > 0 && !is.na(flags) &&
                all(vapply(c("avx2", "fma", "bmi1", "bmi2"),
                           function(f) grepl(paste0("\\b", f, "\\b"), flags), TRUE))
        }
    }

    asset <- if (sysname == "Windows") {
        "inla-windows-x86_64.zip"
    } else if (sysname == "Darwin") {
        if (arm) "inla-macos-arm64-portable.tar.gz" else "inla-macos-x86_64-portable.tar.gz"
    } else if (sysname == "Linux") {
        if (arm) {
            "inla-linux-arm64-armv82-portable.tar.gz"
        } else if (isTRUE(v3) && !is.na(glibc) && glibc >= 2.38) {
            "inla-linux-x86_64-v3-portable.tar.gz"
        } else {
            "inla-linux-x86_64-portable.tar.gz"
        }
    } else {
        stop("no pre-built binary for ", sysname)
    }

    say("platform: ", sysname, " ", machine,
        if (!is.na(glibc)) paste0(", glibc ", glibc) else "",
        if (sysname == "Linux" && !arm) paste0(", x86-64-v3: ", v3) else "")
    say("binary:   ", asset)

    url <- if (is.null(tag)) {
        paste0("https://github.com/", repo, "/releases/latest/download/", asset)
    } else {
        paste0("https://github.com/", repo, "/releases/download/", tag, "/", asset)
    }

    if (is.null(dir)) {
        dir <- file.path(tools::R_user_dir("INLA", "cache"), "stiles-binary",
                         if (is.null(tag)) "latest" else tag)
    }
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    exe <- if (sysname == "Windows") "inla.exe" else "inla"
    bin <- c(Sys.glob(file.path(dir, "*", "bin", exe)), Sys.glob(file.path(dir, "bin", exe)))

    if (length(bin) == 0L || force) {
        arc <- file.path(dir, asset)
        say("downloading ", url)
        utils::download.file(url, arc, mode = "wb", quiet = !verbose)
        say("unpacking")
        if (grepl("\\.zip$", asset)) {
            utils::unzip(arc, exdir = dir)
        } else {
            utils::untar(arc, exdir = dir)
        }
        unlink(arc)
        bin <- c(Sys.glob(file.path(dir, "*", "bin", exe)), Sys.glob(file.path(dir, "bin", exe)))
        if (length(bin) == 0L) stop("no '", exe, "' found under ", dir)
        Sys.chmod(bin[1], "0755")
    } else {
        say("already installed (use force=TRUE to re-download)")
    }

    ## Run it before trusting it: a binary that unpacks but cannot start is the
    ## failure worth catching here, not at the first inla() call.
    ping <- tryCatch(system2(bin[1], "-ping", stdout = TRUE, stderr = TRUE),
                     error = function(e) character(0))
    if (!any(grepl("ALIVE", ping))) {
        warning("the binary did not answer '-ping': ", paste(ping, collapse = " "))
    } else {
        say("binary answers -ping")
    }

    inla.setOption(inla.call = bin[1])
    if (isTRUE(smtp)) inla.setOption(smtp = "stiles")
    say("inla.call = ", bin[1], if (isTRUE(smtp)) ", smtp = stiles" else "")

    info <- file.path(dirname(dirname(bin[1])), "BUILDINFO")
    if (file.exists(info) && verbose) {
        say("BUILDINFO:")
        cat(paste0("    ", readLines(info, warn = FALSE)), sep = "\n")
    }

    invisible(bin[1])
}
