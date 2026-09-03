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
#' version, and the instruction sets the processor actually
#' supports. Each bundle embeds a matching `libstiles`, and carries a
#' `BUILDINFO` file recording the compiler, flags and library versions that
#' produced it.
#'
#' @param tag     Release to install from, e.g. `"Version_26.08.31"` (tags
#'                published before that use the older `"v26.08.18"` form,
#'                and both still resolve). Defaults to the
#'                current release. Pin it for reproducibility.
#' @param dir     Where to unpack. Defaults to a per-user cache directory.
#' @param force   Re-download even when the binary is already installed.
#' @param smtp    Also select the sTiles sparse-matrix backend (default `TRUE`).
#' @param persist Write the settings into `~/.Rprofile` so they survive an R
#'                restart (default `TRUE`). The block is delimited and is
#'                REPLACED on each install, so installing another release
#'                re-points it rather than adding a second entry. A `.bak`
#'                is kept the first time the file is modified.
#' @param repo    GitHub repository holding the releases.
#' @param verbose Report each step.
#'
#' @return The path of the installed binary, invisibly.
#'
#' @examples
#' \dontrun{
#' inla.stiles.install()                  # newest release
#' inla.stiles.install(tag = "Version_26.08.31") # a specific one
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
                                  persist = TRUE,
                                  repo = "hrue/r-inla",
                                  verbose = TRUE) {
    say <- function(...) if (verbose) cat("*", paste0(..., collapse = ""), "\n")

    sysname <- Sys.info()["sysname"]
    machine <- Sys.info()["machine"]
    arm <- machine %in% c("aarch64", "arm64")

    ## The C library decides which Linux bundle can run at all, and the
    ## processor decides which one is safe: x86-64-v3 needs AVX2, FMA and BMI,
    ## meaning Haswell (2013) or newer. A recent distribution on an older
    ## processor is an ordinary combination, and choosing by glibc alone hands
    ## it a binary that aborts on the first vectorised kernel.
    ## Every bundle has an instruction-set floor, so a machine below it gets a
    ## binary that aborts with an illegal instruction on the first vectorised
    ## kernel. Read the floors here and refuse rather than install something
    ## that cannot run.
    cpuflags <- function(key) {
        if (!file.exists("/proc/cpuinfo")) return(NA_character_)
        ln <- grep(paste0("^", key), readLines("/proc/cpuinfo", warn = FALSE), value = TRUE)
        if (!length(ln)) NA_character_ else ln[1]
    }
    has <- function(line, want) {
        !is.na(line) && all(vapply(want,
            function(f) grepl(paste0("\\b", f, "\\b"), line), TRUE))
    }

    glibc <- NA_real_
    v3 <- FALSE
    avx2 <- NA          # NA when it cannot be determined; do not block on that
    armv82 <- NA
    if (sysname == "Linux") {
        s <- tryCatch(system("getconf GNU_LIBC_VERSION", intern = TRUE),
                      error = function(e) NA_character_,
                      warning = function(w) NA_character_)
        if (!is.na(s[1])) glibc <- suppressWarnings(as.numeric(strsplit(s, " ")[[1]][2]))
        if (arm) {
            ## armv8.2 markers: the half-precision and dot-product extensions
            ## the arm64 bundle is compiled for. Present on Graviton2 and newer,
            ## Ampere, Grace and the Raspberry Pi 5; absent on armv8.0 parts such
            ## as the Raspberry Pi 4 and Graviton1.
            f <- cpuflags("Features")
            if (!is.na(f)) armv82 <- has(f, c("asimdrdm")) && has(f, c("fphp", "asimdhp"))
        } else {
            f <- cpuflags("flags")
            if (!is.na(f)) {
                avx2 <- has(f, "avx2")
                v3 <- has(f, c("avx2", "fma", "bmi1", "bmi2"))
            }
        }
    }

    asset <- if (sysname == "Windows") {
        "inla-windows-x86_64.zip"
    } else if (sysname == "Darwin") {
        ## The macOS bundles carry libraries built on a recent system, so an
        ## older one can fail at load time with a message that explains nothing.
        mac <- tryCatch(system("sw_vers -productVersion", intern = TRUE)[1],
                        error = function(e) NA_character_,
                        warning = function(w) NA_character_)
        if (!is.na(mac)) {
            major <- suppressWarnings(as.numeric(strsplit(mac, "\\.")[[1]][1]))
            if (!is.na(major) && major < 11) {
                warning("macOS ", mac, " is older than the bundles are built for ",
                        "(11.0). It may fail to load; see the release page.")
            }
        }
        if (arm) "inla-macos-arm64-portable.tar.gz" else "inla-macos-x86_64-portable.tar.gz"
    } else if (sysname == "Linux") {
        if (arm) {
            if (identical(armv82, FALSE)) {
                stop("this machine is armv8.0 (a Raspberry Pi 4 or Graviton1, say). ",
                     "The only arm64 build needs armv8.2, so no pre-built binary ",
                     "fits it; build from source instead.")
            }
            "inla-linux-arm64-armv82-portable.tar.gz"
        } else if (isTRUE(v3) && !is.na(glibc) && glibc >= 2.38) {
            "inla-linux-x86_64-v3-portable.tar.gz"
        } else {
            if (identical(avx2, FALSE)) {
                stop("this processor has no AVX2, which every x86-64 build needs ",
                     "(Haswell 2013 or newer). Build from source instead.")
            }
            "inla-linux-x86_64-portable.tar.gz"
        }
    } else {
        stop("no pre-built binary for ", sysname)
    }

    say("platform: ", sysname, " ", machine,
        if (!is.na(glibc)) paste0(", glibc ", glibc) else "",
        if (sysname == "Linux" && !arm) paste0(", x86-64-v3: ", v3) else "",
        if (sysname == "Linux" && arm && !is.na(armv82)) paste0(", armv8.2: ", armv82) else "")
    say("binary:   ", asset)

    ## Resolve which release "latest" currently means, and cache under THAT
    ## name. Caching under a directory called "latest" makes the first install
    ## permanent: a newer release lands on the same path, the binary is found,
    ## and the caller keeps the old one until they think to pass force = TRUE.
    ## One request to the releases API, and the tag comes back in the JSON; if
    ## it cannot be reached, fall back to the redirecting URL as before.
    resolved <- tag
    if (is.null(resolved)) {
        js <- tryCatch(
            paste(readLines(paste0("https://api.github.com/repos/", repo,
                                   "/releases/latest"), warn = FALSE),
                  collapse = ""),
            error = function(e) "", warning = function(w) "")
        m <- regmatches(js, regexpr('"tag_name"[^"]*"[^"]+"', js))
        if (length(m) == 1L) {
            resolved <- sub('.*"tag_name"[^"]*"([^"]+)".*', "\\1", m)
            say("release:  ", resolved)
        }
    }

    url <- if (is.null(tag) && is.null(resolved)) {
        paste0("https://github.com/", repo, "/releases/latest/download/", asset)
    } else {
        paste0("https://github.com/", repo, "/releases/download/",
               if (is.null(tag)) resolved else tag, "/", asset)
    }

    default.dir <- is.null(dir)
    if (is.null(dir)) {
        dir <- file.path(tools::R_user_dir("INLA", "cache"), "stiles-binary",
                         if (!is.null(resolved)) resolved else "latest")
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

    ## The bundles ship bin/inla.run beside the binary (it preloads the
    ## bundled allocator, then execs inla); that script is the upstream entry
    ## point, so point inla.call at it when it exists. Older releases have
    ## only the binary, and Windows has no wrapper.
    if (sysname != "Windows") {
        run <- file.path(dirname(bin[1]), "inla.run")
        if (file.exists(run)) {
            Sys.chmod(run, "0755")
            bin[1] <- run
        }
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

    ## A stable "latest" entry beside the versioned ones (mirrors ~/R/*/default
    ## and MKL's own .../mkl/latest), so a path saved in ~/.Rprofile survives a
    ## new release instead of needing hand-editing every time one lands. Only
    ## when following the latest release (a pinned tag means the caller wants
    ## an exact, unmoving path) -- and never when `dir` itself is ALREADY named
    ## "latest", which happens when the release-tag API call failed but the
    ## download redirect still worked (see `resolved` above): there is nothing
    ## to point at that is not already sitting at that name.
    alive <- any(grepl("ALIVE", ping))
    if (default.dir && is.null(tag) && alive && !identical(basename(dir), "latest")) {
        latest <- file.path(dirname(dir), "latest")
        if (file.exists(latest) || nzchar(Sys.readlink(latest)))
            unlink(latest, recursive = TRUE, force = TRUE)
        ## Relative target: the link (and the cache root, if the user ever
        ## relocates it) keeps working without repointing.
        made <- tryCatch(file.symlink(basename(dir), latest), error = function(e) FALSE)
        if (!isTRUE(made)) {
            ## Symlinks need a privilege Windows does not always grant; a
            ## real copy costs disk (one release, ~100 MB) but always works.
            ## dir and latest are always siblings (same parent), so copying
            ## `dir` itself INTO dirname(latest) would try to create dir at
            ## its own existing path -- a self-copy that corrupts it. Copy
            ## CONTENTS into a freshly made "latest" instead.
            say("symlink unavailable, copying to 'latest' instead")
            dir.create(latest, recursive = TRUE, showWarnings = FALSE)
            ok <- file.copy(list.files(dir, full.names = TRUE), latest,
                            recursive = TRUE, copy.mode = TRUE)
            if (!all(ok)) warning("copying to 'latest' was incomplete")
        }
        ## bin[1] already reflects the inla.run substitution above; keep that
        ## choice, just reached through the stable name instead of the
        ## version-pinned one.
        bin[1] <- file.path(latest, substring(bin[1], nchar(dir) + 2L))
        say("stable path: ", bin[1])
    }

    ## Superseded releases serve nobody once a newer one has answered -ping:
    ## each is ~100 MB, and keeping them is how a cache quietly grows to
    ## gigabytes. Clean AFTER the ping, never before, and only in the default
    ## cache (a caller-supplied dir has siblings that are none of our
    ## business) when following the latest release (a pinned tag means the
    ## caller manages versions deliberately). "latest" is never superseded --
    ## it was just re-pointed above, or (the offline-fallback case) it IS the
    ## live install -- either way it must survive this pass.
    if (default.dir && is.null(tag) && alive) {
        for (d in list.dirs(dirname(dir), recursive = FALSE)) {
            if (!identical(basename(d), basename(dir)) && !identical(basename(d), "latest")) {
                say("removing superseded ", basename(d))
                unlink(d, recursive = TRUE)
            }
        }
    }

    inla.setOption(inla.call = bin[1])
    if (isTRUE(smtp)) inla.setOption(smtp = "stiles")
    say("inla.call = ", bin[1], if (isTRUE(smtp)) ", smtp = stiles" else "")

    info <- file.path(dirname(dirname(bin[1])), "BUILDINFO")
    if (file.exists(info) && verbose) {
        say("BUILDINFO:")
        cat(paste0("    ", readLines(info, warn = FALSE)), sep = "\n")
    }

    ## The options set above live in this R session only; a restart keeps the
    ## downloaded binary but forgets both settings, which looks like the
    ## install vanished. Print the calls that bring it back, unconditionally:
    ## they are the actionable output of this function, not progress chatter.
    if (isTRUE(persist)) {
        ## Write the settings into ~/.Rprofile so a restart keeps them.
        ##
        ## The block is delimited by markers and REPLACED on every install, so
        ## installing a second release re-points the same lines instead of
        ## appending a second inla.call that would shadow the first depending
        ## on order. Everything outside the markers is copied through
        ## untouched, and a one-time .bak is kept the first time this file is
        ## modified, because it is the user's file and not ours to lose.
        ok <- tryCatch({
            rp <- path.expand("~/.Rprofile")
            beg <- "## >>> INLA: set by inla.stiles.install() >>>"
            end <- "## <<< INLA <<<"
            body <- c(beg,
                      sprintf('INLA::inla.setOption(inla.call = "%s")', bin[1]),
                      if (isTRUE(smtp)) 'INLA::inla.setOption(smtp = "stiles")',
                      end)
            oldl <- if (file.exists(rp)) readLines(rp, warn = FALSE) else character(0)
            if (length(oldl) && !file.exists(paste0(rp, ".bak"))) {
                try(writeLines(oldl, paste0(rp, ".bak")), silent = TRUE)
            }
            i <- which(oldl == beg); j <- which(oldl == end)
            keep <- if (length(i) == 1L && length(j) == 1L && j > i) {
                c(utils::head(oldl, i - 1L), utils::tail(oldl, length(oldl) - j))
            } else {
                oldl
            }
            ## Write beside and rename: an interrupted write must not leave a
            ## truncated .Rprofile, which would break every future R session.
            tmp <- paste0(rp, ".new")
            writeLines(c(keep, body), tmp)
            file.rename(tmp, rp)
            TRUE
        }, error = function(e) e)

        if (isTRUE(ok)) {
            say("~/.Rprofile updated; the settings survive a restart")
            cat("\nWritten to ~/.Rprofile (replacing any previous INLA block):\n\n")
            cat(sprintf('    INLA::inla.setOption(inla.call = "%s")\n', bin[1]))
            if (isTRUE(smtp)) cat('    INLA::inla.setOption(smtp = "stiles")\n')
            cat("\nRun with persist = FALSE to leave ~/.Rprofile alone.\n")
        } else {
            cat("\nCould not update ~/.Rprofile (", conditionMessage(ok), ").\n", sep = "")
            cat("Add these lines yourself to make the settings permanent:\n\n")
            cat(sprintf('    INLA::inla.setOption(inla.call = "%s")\n', bin[1]))
            if (isTRUE(smtp)) cat('    INLA::inla.setOption(smtp = "stiles")\n')
        }
    } else {
        cat("\nThese settings do not survive an R restart. To restore them next time, run:\n\n")
        cat(sprintf('    inla.setOption(inla.call = "%s")\n', bin[1]))
        if (isTRUE(smtp)) cat('    inla.setOption(smtp = "stiles")\n')
        cat("\nor re-run with persist = TRUE to write them to ~/.Rprofile.\n")
    }

    invisible(bin[1])
}


#' @title What the session is using right now
#'
#' @description
#' `inla.stiles.status()` answers the question `inla.stiles.install()` leaves
#' open after a restart: which binary is this session actually pointing at,
#' with which backend, and what was that binary built from. Everything is read
#' from the session's options and the binary's own `BUILDINFO`, not from what
#' an earlier install intended.
#'
#' @param ping Also run the binary with `-ping` to prove it starts (default
#'             `TRUE`; costs about a second).
#' @param verbose Print the report (default `TRUE`).
#'
#' @return Invisibly, a list with `inla.call`, `smtp`, `release`, `alive`,
#'         `buildinfo` (character vector) and `cache` (installed releases).
#'
#' @seealso [inla.stiles.install()]
#' @rdname stiles.install
#' @export inla.stiles.status

`inla.stiles.status` <- function(ping = TRUE, verbose = TRUE) {
    say <- function(...) if (verbose) cat("*", paste0(..., collapse = ""), "\n")

    call <- tryCatch(inla.getOption("inla.call"), error = function(e) NULL)
    smtp <- tryCatch(inla.getOption("smtp"), error = function(e) NULL)

    cache <- file.path(tools::R_user_dir("INLA", "cache"), "stiles-binary")
    installed <- basename(list.dirs(cache, recursive = FALSE))

    ## Is the active binary one of ours, and which release? The cache path is
    ## <cache>/<tag>/bin/<exe>, so the tag is two levels up from the file.
    release <- NA_character_
    if (!is.null(call) && is.character(call) && nzchar(call) && file.exists(call)) {
        root <- dirname(dirname(call))
        if (identical(normalizePath(dirname(root), mustWork = FALSE),
                      normalizePath(cache, mustWork = FALSE))) {
            release <- basename(root)
        }
        info <- file.path(root, "BUILDINFO")
        buildinfo <- if (file.exists(info)) readLines(info, warn = FALSE) else character(0)
    } else {
        buildinfo <- character(0)
    }

    say("inla.call: ", if (is.null(call) || !is.character(call) || !nzchar(call[1])) "(package default)" else call)
    if (!is.na(release)) say("release:   ", release, "  (from the install cache)")
    say("smtp:      ", if (is.null(smtp)) "(default)" else smtp)

    alive <- NA
    if (isTRUE(ping) && !is.null(call) && is.character(call) && file.exists(call)) {
        out <- tryCatch(system2(call, "-ping", stdout = TRUE, stderr = TRUE),
                        error = function(e) character(0))
        alive <- any(grepl("ALIVE", out))
        say(if (isTRUE(alive)) "binary answers -ping" else "binary DID NOT answer -ping")
    }

    ## The sTiles release this binary carries. BUILDINFO records it as
    ## "libstiles:  Version_2026.09.03"; that tag is what the sTiles release
    ## page and the pairing with an INLA bundle are keyed on, so report it as
    ## a field rather than leaving the caller to parse buildinfo themselves.
    stiles.version <- NA_character_
    if (length(buildinfo)) {
        ln <- grep("^[[:space:]]*libstiles[[:space:]]*:", buildinfo,
                   ignore.case = TRUE, value = TRUE)
        if (length(ln)) {
            v <- sub("^[[:space:]]*[Ll]ibstiles[[:space:]]*:[[:space:]]*", "", ln[1])
            v <- trimws(v)
            if (nzchar(v)) stiles.version <- v
        }
    }
    say("sTiles:    ", if (is.na(stiles.version)) "(unknown; no BUILDINFO)" else stiles.version)
    if (length(buildinfo) && verbose) {
        keep <- grep("^(compiler|libstiles|blas|date):", buildinfo, ignore.case = TRUE, value = TRUE)
        if (length(keep)) { say("BUILDINFO:"); cat(paste0("    ", keep), sep = "\n") }
    }
    if (length(installed)) say("installed releases: ", paste(installed, collapse = ", "))

    invisible(list(inla.call = call, smtp = smtp, release = release,
                   stiles.version = stiles.version,
                   alive = alive, buildinfo = buildinfo, cache = installed))
}
