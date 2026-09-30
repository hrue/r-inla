## Can we actually write there? Tested by writing, not by file.access(): on a
## network filesystem or under an ACL the permission bits say yes and the write
## still fails, which would send the install to a directory it cannot fill.
## Cleans up after itself and leaves no empty tree behind on failure.
`inla.dir.writable` <- function(dir) {
    made <- !dir.exists(dir)
    if (made && !dir.create(dir, recursive = TRUE, showWarnings = FALSE)) {
        return(FALSE)
    }
    probe <- file.path(dir, paste0(".inla-write-test-", Sys.getpid()))
    ok <- isTRUE(tryCatch({
        cat("", file = probe)
        file.exists(probe)
    }, error = function(e) FALSE, warning = function(w) FALSE))
    unlink(probe)
    if (made && !ok) unlink(dir, recursive = TRUE)
    ok
}

## The repository holding the releases. One value, overridable for a fork or
## a move, so the name is not spelled out in five places.
`inla.repo` <- function() Sys.getenv("INLA_REPO", "hrue/r-inla")

## Per-user cache directory for downloaded binaries.
##
## tools::R_user_dir() arrived in R 4.0 but this package supports R >= 3.5, so
## calling it directly would error on exactly the older installations most
## likely to want a prebuilt binary. Use it when present, and otherwise the
## same conventional location it would return on unix.
`inla.cache.dir` <- function() {
    if (exists("R_user_dir", envir = asNamespace("tools"), inherits = FALSE)) {
        return(get("R_user_dir", envir = asNamespace("tools"))("INLA", "cache"))
    }
    base <- Sys.getenv("XDG_CACHE_HOME")
    if (!nzchar(base)) {
        base <- if (.Platform$OS.type == "windows") {
            file.path(Sys.getenv("LOCALAPPDATA"), "R", "cache")
        } else {
            path.expand("~/.cache")
        }
    }
    file.path(base, "R", "INLA")
}

#' @title Install a pre-built INLA binary with sTiles support
#'
#' @description
#' `inla.stiles.install()` works out which pre-built binary matches this
#' machine, downloads the release this package version needs, unpacks it
#' beside the package, and checks that it runs. Nothing to configure
#' afterwards: `inla.call` stays unset and the binary is found where it was
#' put.
#'
#' The binaries are portable: six of them cover Linux (x86-64 baseline and
#' x86-64-v3), Linux arm64, macOS (Intel and Apple silicon) and Windows, so the
#' choice follows from the operating system, the architecture, the C library
#' version, and the instruction sets the processor actually
#' supports. Each bundle embeds a matching `libstiles`, and carries a
#' `BUILDINFO` file recording the compiler, flags and library versions that
#' produced it.
#'
#' @param dir     Where to unpack. By default the installed package, so the
#'                binary travels with it, or a per-user cache when the package
#'                directory is not writable. Those two are the only places
#'                INLA searches, so a binary put anywhere else has to be named
#'                with `inla.setOption(inla.call = ...)` before anything uses
#'                it; the function prints that line.
#' @param force   Re-download even when the binary is already installed.
#' @param smtp    Also select the sTiles sparse-matrix backend (default `TRUE`).
#' @param verbose Report each step.
#'
#' @return The path of the installed binary, invisibly.
#'
#' @examples
#' \dontrun{
#' inla.stiles.install()   # the binary this package version needs
#' }
#'
#' @seealso [inla.binary.install()]
#' @name stiles.install
#' @aliases inla.stiles.install
#' @rdname stiles.install
#' @export inla.stiles.install

`inla.stiles.install` <- function(dir = NULL,
                                  force = FALSE,
                                  smtp = TRUE,
                                  verbose = TRUE) {
    say <- function(...) if (verbose) cat("*", paste0(..., collapse = ""), "\n")
    repo <- inla.repo()

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
    tag_name_of <- function(u) {
        js <- tryCatch(paste(readLines(u, warn = FALSE), collapse = ""),
                       error = function(e) "", warning = function(w) "")
        m <- regmatches(js, regexpr('"tag_name"[^"]*"[^"]+"', js))
        if (length(m) == 1L) sub('.*"tag_name"[^"]*"([^"]+)".*', "\\1", m) else NULL
    }

    resolved <- NULL
    {
        ## Default to the release that MATCHES THIS R PACKAGE, not merely the
        ## newest one. The binary and the package are two halves of one
        ## release: pairing 26.09.03 with a binary from a later release is how
        ## a user ends up debugging a mismatch that nobody intended. The
        ## release tag is the package's Version with a "v" in front, read from
        ## DESCRIPTION rather than packageVersion(), because R normalises
        ## "26.09.03" to "26.9.3" and the tag keeps the zeros.
        ## Which BINARY release this R package needs, declared in DESCRIPTION
        ## as Config/INLA/BinaryVersion. It is now always the SAME string as
        ## the package's own Version: one number identifies the R package and
        ## the binary that belongs with it.
        ##
        ## It used to move only when the C sources changed, so an R-only fix
        ## did not force a binary release. That left two similar-looking dates
        ## that disagreed (Version 26.09.07-1 against BinaryVersion 26.09.07),
        ## and no way to say which one identified what a user had.
        ##
        ## The rule this creates: every release must publish binaries, since
        ## this field names a release tag that has to exist.
        ##
        ## Read through this field rather than Version anyway, so an
        ## installation predating the change still works, and so the pairing
        ## has one authority.
        pv <- tryCatch(utils::packageDescription("INLA")[["Config/INLA/BinaryVersion"]],
                       error = function(e) NULL)
        if (is.null(pv) || !nzchar(pv)) {
            pv <- tryCatch(utils::packageDescription("INLA")$Version,
                           error = function(e) NULL)
        }
        if (!is.null(pv) && nzchar(pv)) {
            for (cand in c(paste0("v", pv), paste0("Version_", pv))) {
                if (!is.null(tag_name_of(paste0("https://api.github.com/repos/", repo,
                                                "/releases/tags/", cand)))) {
                    resolved <- cand
                    say("release:  ", resolved, "  (matches this R package)")
                    break
                }
            }
        }
        ## No release for this package version (a development build, or one
        ## whose release was never cut): fall back to the newest, and say so,
        ## because that pairing is then not guaranteed.
        if (is.null(resolved)) {
            resolved <- tag_name_of(paste0("https://api.github.com/repos/", repo,
                                           "/releases/latest"))
            if (!is.null(resolved)) {
                say("release:  ", resolved,
                    "  (NO release matches package ", if (is.null(pv)) "?" else pv,
                    "; using the newest)")
            }
        }
    }

    url <- if (is.null(resolved)) {
        paste0("https://github.com/", repo, "/releases/latest/download/", asset)
    } else {
        paste0("https://github.com/", repo, "/releases/download/", resolved, "/", asset)
    }

    default.dir <- is.null(dir)
    if (is.null(dir)) {
        ## Prefer the installed package: the binary then travels with it, so a
        ## NULL inla.call finds the one that belongs to this version and no
        ## cache has to be kept in step. The cache is the fallback for a site
        ## library, where the package directory belongs to root. Reinstalling
        ## the package removes a package-side binary, which is accepted.
        tag.dir <- if (!is.null(resolved)) resolved else "latest"
        pkg <- tryCatch(find.package("INLA"), error = function(e) "")
        cand <- if (nzchar(pkg)) file.path(pkg, "stiles-binary", tag.dir) else ""
        dir <- ""
        if (nzchar(cand) && inla.dir.writable(cand)) dir <- cand
        if (!nzchar(dir)) {
            dir <- file.path(inla.cache.dir(), "stiles-binary", tag.dir)
        }
        say("target:   ", dir,
            if (identical(dir, cand)) "  (inside the package)" else "  (user cache)")
    }
    dir.create(dir, recursive = TRUE, showWarnings = FALSE)

    exe <- if (sysname == "Windows") "inla.exe" else "inla"

    ## Where the binary lands depends on how the bundle was packed, and the
    ## layouts differ per platform. The unix bundles use <root>/bin/inla, and
    ## unpack either into dir/ or into a single top directory inside it. The
    ## WINDOWS zip is flat: inla.exe and its DLLs sit directly in dir, with no
    ## bin/ at all, so the two bin/ patterns alone found nothing and the
    ## install failed with "no 'inla.exe' found" while the file was plainly
    ## there. Search all four shapes, nearest first.
    find_bin <- function(d) {
        c(Sys.glob(file.path(d, "*", "bin", exe)),
          Sys.glob(file.path(d, "bin", exe)),
          Sys.glob(file.path(d, "*", exe)),
          Sys.glob(file.path(d, exe)))
    }
    bin <- find_bin(dir)

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
        bin <- find_bin(dir)
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
    if (default.dir && alive && !identical(basename(dir), "latest")) {
        latest <- file.path(dirname(dir), "latest")
        ## Remove whatever is there, symlink or directory, and VERIFY it went.
        ##
        ## On Windows a directory symlink is a reparse point, and R's own
        ## unlink() can fail on it:
        ##     cannot delete reparse point ... mismatch between the tag
        ##     specified in the request and the tag present
        ## It warns and returns as if nothing happened, so the stale link
        ## survived and the copy below wrote into a symlink, dying with the
        ## unhelpful "more 'from' files than 'to' files".
        ##
        ## Which call succeeds depends on how the link was made (symlink or
        ## junction) and on the Windows version, so try the plausible ones in
        ## order and CHECK after each rather than trusting any single one.
        ## Nothing here can be tested on Linux: the failure needs a reparse
        ## point, so the Windows lane is the only real test.
        ## Sys.readlink() returns "" for a real file, the target for a link,
        ## and NA when the path does not exist. nzchar(NA) is TRUE, so the
        ## obvious !nzchar(Sys.readlink(x)) test reports "still there" for a
        ## path that is already gone. Handle the NA explicitly.
        gone <- function() {
            rl <- suppressWarnings(Sys.readlink(latest))
            is_link <- !is.na(rl) && nzchar(rl)
            !file.exists(latest) && !is_link
        }
        if (!gone()) suppressWarnings(unlink(latest, recursive = TRUE, force = TRUE))
        if (!gone()) suppressWarnings(unlink(latest, force = TRUE))
        if (!gone()) suppressWarnings(try(file.remove(latest), silent = TRUE))
        if (!gone()) {
            stop("could not replace the stable path '", latest, "'.\n",
                 "  Delete it by hand and re-run: unlink(\"", latest, "\")\n",
                 "  or, on Windows: rmdir \"", gsub("/", "\\\\", latest), "\"")
        }
        ## Relative target: the link (and the cache root, if the user ever
        ## relocates it) keeps working without repointing.
        ## suppressWarnings, not just tryCatch: on Windows file.symlink()
        ## does not raise an error, it WARNS and returns FALSE. The warning is
        ## deferred to the end of the session, so a user saw
        ##     cannot symlink 'Version_...' to '.../latest', reason
        ##     'A required privilege is not held by the client'
        ## printed after a successful install, which reads like a failure. The
        ## copy below handles it; there is nothing for anyone to act on.
        ## INLA_STILES_NO_SYMLINK forces the copy branch below. It exists for
        ## CI: the GitHub Windows runner is an ADMINISTRATOR (runneradmin) and
        ## therefore holds SeCreateSymbolicLinkPrivilege, so the symlink always
        ## succeeds there and the fallback that every ordinary Windows user
        ## takes was never once executed in a test.
        made <- if (nzchar(Sys.getenv("INLA_STILES_NO_SYMLINK"))) {
            FALSE
        } else {
            suppressWarnings(
                tryCatch(file.symlink(basename(dir), latest), error = function(e) FALSE))
        }
        if (!isTRUE(made)) {
            ## Symlinks need a privilege Windows does not always grant; a
            ## real copy costs disk (one release, ~100 MB) but always works.
            ## dir and latest are always siblings (same parent), so copying
            ## `dir` itself INTO dirname(latest) would try to create dir at
            ## its own existing path -- a self-copy that corrupts it. Copy
            ## CONTENTS into a freshly made "latest" instead.
            say("symlink unavailable, copying to 'latest' instead")
            dir.create(latest, recursive = TRUE, showWarnings = FALSE)
            ## file.copy() needs `to` to be one EXISTING directory, otherwise it
            ## pairs from/to elementwise and fails with "more 'from' files than
            ## 'to' files", which says nothing about the real cause. Check it.
            if (!dir.exists(latest)) {
                stop("could not create the stable path '", latest, "'")
            }
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
    if (default.dir && alive) {
        for (d in list.dirs(dirname(dir), recursive = FALSE)) {
            if (!identical(basename(d), basename(dir)) && !identical(basename(d), "latest")) {
                say("removing superseded ", basename(d))
                unlink(d, recursive = TRUE)
            }
        }
    }

    ## inla.call is left NULL on purpose. NULL means "the binary that belongs to
    ## this package version", which the resolver finds at the paths written
    ## above and which stays right after an upgrade. Storing the path here made
    ## every install an override pinned to one release. Setting it is the
    ## expert's choice, not the installer's.
    if (isTRUE(smtp)) inla.setOption(smtp = "stiles")
    say("installed: ", bin[1], if (isTRUE(smtp)) "   smtp = stiles" else "")

    ## BUILDINFO sits at the bundle root, which is the binary's grandparent for
    ## a <root>/bin/inla layout but the binary's OWN directory for the flat
    ## Windows zip. Take whichever exists rather than assuming the depth.
    info <- c(file.path(dirname(dirname(bin[1])), "BUILDINFO"),
              file.path(dirname(bin[1]), "BUILDINFO"))
    info <- c(info[file.exists(info)], info[1])[1]
    if (file.exists(info) && verbose) {
        say("BUILDINFO:")
        cat(paste0("    ", readLines(info, warn = FALSE)), sep = "\n")
    }

    ## Nothing to remember: inla.call stays unset, so the binary just
    ## installed is the one found from now on, here and after a restart.
    ##
    ## Unless `dir` sent it somewhere else. The resolver only searches the
    ## package and the cache, so a binary anywhere else is installed and then
    ## ignored, silently and especially after a restart. Say so, and give the
    ## line that makes it usable, rather than letting it look like it worked.
    if (!default.dir) {
        cat("\nThis is outside the places INLA searches, so nothing will use it\n")
        cat("until you say so:\n\n")
        cat(sprintf('    inla.setOption(inla.call = "%s")\n', bin[1]))
        cat("\nPut that line in ~/.Rprofile to keep it across restarts, or run\n")
        cat("inla.stiles.install() with no dir to install where INLA looks.\n")
    } else {
        cat("\nDone. Check with inla.stiles.status().\n")
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

    ## RESOLVE, do not read the option. inla.call is NULL on a healthy
    ## installation, so reading it reported "(package default)" and then said
    ## nothing else: no release, no ping, no BUILDINFO, because every field
    ## below is keyed on this path. Ask for the binary that will actually run.
    set <- tryCatch(inla.getOption("inla.call"), error = function(e) NULL)
    set <- if (is.null(set) || !is.character(set) || !nzchar(set[1])) NULL else set[1]
    call <- tryCatch(inla.binary.path(check = FALSE), error = function(e) NA_character_)
    if (length(call) != 1L || is.na(call) || !nzchar(call)) call <- NULL
    smtp <- tryCatch(inla.getOption("smtp"), error = function(e) NULL)

    ## Both places a binary can live, and everything installed in either.
    pkg <- tryCatch(find.package("INLA"), error = function(e) "")
    pkg.base <- if (nzchar(pkg)) file.path(pkg, "stiles-binary") else ""
    cache <- file.path(inla.cache.dir(), "stiles-binary")
    installed <- c(
        if (nzchar(pkg.base)) paste0(basename(list.dirs(pkg.base, recursive = FALSE)),
                                     " (package)") else character(0),
        paste0(basename(list.dirs(cache, recursive = FALSE)), " (cache)"))

    ## Which of the two it came from, and which release. Both layouts are
    ## <base>/<tag>/bin/<exe>, so the tag is two levels up and the base three.
    release <- NA_character_
    source <- if (!is.null(set)) "set in inla.call" else NA_character_
    if (!is.null(call) && file.exists(call)) {
        root <- dirname(dirname(call))
        base <- normalizePath(dirname(root), mustWork = FALSE)
        if (identical(base, normalizePath(cache, mustWork = FALSE))) {
            release <- basename(root)
            if (is.na(source)) source <- "user cache"
        } else if (nzchar(pkg.base) &&
                   identical(base, normalizePath(pkg.base, mustWork = FALSE))) {
            release <- basename(root)
            if (is.na(source)) source <- "inside the package"
        } else if (is.na(source)) {
            source <- "elsewhere"
        }
        info <- file.path(root, "BUILDINFO")
        buildinfo <- if (file.exists(info)) readLines(info, warn = FALSE) else character(0)
    } else {
        buildinfo <- character(0)
    }

    say("inla.call: ", if (is.null(call)) "(no binary found)" else call)
    if (!is.na(release)) say("release:   ", release)
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

    ## The four versions that decide whether this installation is coherent,
    ## in one place. They are allowed to differ: the R package moves on its
    ## own, and only a binary OLDER than the declared requirement is a
    ## problem. Printing all four is what makes a mismatch obvious instead of
    ## something a user has to reconstruct from three separate commands.
    r.version <- tryCatch(utils::packageDescription("INLA")$Version,
                          error = function(e) NA_character_)
    binary.required <- tryCatch(utils::packageDescription("INLA")[["Config/INLA/BinaryVersion"]],
                                error = function(e) NULL)
    if (is.null(binary.required) || !nzchar(binary.required)) binary.required <- NA_character_
    binary.version <- NA_character_
    if (!is.null(call) && is.character(call) && nzchar(call) && file.exists(call)) {
        vout <- suppressWarnings(tryCatch(
            system2(call, "-V", stdout = TRUE, stderr = TRUE, timeout = 20),
            error = function(e) character(0)))
        hit <- grep("version", vout, ignore.case = TRUE, value = TRUE)
        if (length(hit)) binary.version <- trimws(sub(".*version:[[:space:]]*", "", hit[1]))
    }
    ## Exact, not "or newer": one number identifies the R package and the
    ## binary that belongs with it, so anything else is a pairing nobody
    ## tested. A binary chosen by hand is the user's business and says so.
    ok <- NA
    if (!is.na(binary.version) && !is.na(binary.required)) {
        ok <- tryCatch(package_version(binary.version) == package_version(binary.required),
                       error = function(e) NA)
    }
    say("R package: ", if (is.na(r.version)) "?" else r.version)
    say("binary:    ", if (is.na(binary.version)) "(none)" else binary.version,
        "   (this package needs ",
        if (is.na(binary.required)) "?" else binary.required, ": ",
        if (is.na(ok)) "unknown"
        else if (isTRUE(ok)) "OK"
        else if (!is.null(set)) "chosen by you"
        else "MISMATCH, run inla.stiles.install()", ")")
    if (length(buildinfo) && verbose) {
        keep <- grep("^(compiler|libstiles|blas|date):", buildinfo, ignore.case = TRUE, value = TRUE)
        if (length(keep)) { say("BUILDINFO:"); cat(paste0("    ", keep), sep = "\n") }
    }
    if (length(installed)) say("installed releases: ", paste(installed, collapse = ", "))

    invisible(list(inla.call = call, smtp = smtp, release = release,
                   stiles.version = stiles.version,
                   r.version = r.version,
                   binary.version = binary.version,
                   binary.required = binary.required,
                   binary.ok = ok,
                   alive = alive, buildinfo = buildinfo, cache = installed))
}

#' @title Binary releases available to install
#'
#' @description
#' List the `inla` binary releases published for [inla.stiles.install()],
#' marking which are already on this machine, which one is in use, and which
#' one this R package asks for. The companion to [inla.stiles.status()], which
#' reports only what is installed.
#'
#' @param n     How many of the most recent releases to list.
#' @param verbose Print the table.
#'
#' @returns Invisibly, a `data.frame` with `tag`, `published`, `installed`,
#'   `active` and `required`.
#'
#' @examples
#' \dontrun{
#' inla.stiles.releases()
#' subset(inla.stiles.releases(verbose = FALSE), installed)
#' }
#'
#' @seealso [inla.stiles.install()], [inla.stiles.status()]
#' @export inla.stiles.releases

`inla.stiles.releases` <- function(n = 10, verbose = TRUE) {
    repo <- inla.repo()
    ## Same dependency-free parse as the installer: three flat fields out of
    ## the releases JSON, rather than adding jsonlite to Imports for this.
    js <- tryCatch(paste(readLines(paste0("https://api.github.com/repos/", repo,
                                          "/releases?per_page=", as.integer(n)),
                                   warn = FALSE), collapse = " "),
                   error = function(e) "", warning = function(w) "")
    if (!nzchar(js)) {
        if (verbose) cat("  Could not reach the releases index for ", repo, "\n", sep = "")
        return(invisible(data.frame()))
    }
    tags <- regmatches(js, gregexpr('"tag_name"[[:space:]]*:[[:space:]]*"[^"]+"', js))[[1]]
    tags <- sub('.*"([^"]+)"$', "\\1", tags)
    pub  <- regmatches(js, gregexpr('"published_at"[[:space:]]*:[[:space:]]*"[^"]+"', js))[[1]]
    pub  <- substr(sub('.*"([^"]+)"$', "\\1", pub), 1, 10)
    if (length(pub) != length(tags)) pub <- rep(NA_character_, length(tags))
    if (!length(tags)) return(invisible(data.frame()))

    ## Installed in EITHER place, not just the cache: the default target is
    ## now the package itself.
    pkg <- tryCatch(find.package("INLA"), error = function(e) "")
    have <- basename(c(
        if (nzchar(pkg)) list.dirs(file.path(pkg, "stiles-binary"), recursive = FALSE)
        else character(0),
        list.dirs(file.path(inla.cache.dir(), "stiles-binary"), recursive = FALSE)))
    ## Resolve, do not read the option: it is NULL on a healthy installation,
    ## so reading it marked nothing as in use.
    call <- tryCatch(inla.binary.path(check = FALSE), error = function(e) NULL)
    ## The active release is the directory the running binary sits in:
    ## <base>/<tag>/bin/<exe>, so the tag is two levels up.
    active <- NA_character_
    if (length(call) == 1L && !is.na(call) && nzchar(call)) {
        active <- basename(dirname(dirname(call)))
    }
    need <- tryCatch(utils::packageDescription("INLA")[["Config/INLA/BinaryVersion"]],
                     error = function(e) NULL)
    if (is.null(need)) need <- NA_character_

    out <- data.frame(tag = tags, published = pub,
                      installed = tags %in% have,
                      active = !is.na(active) & tags == active,
                      required = !is.na(need) &
                          (tags == paste0("v", need) | tags == paste0("Version_", need)),
                      stringsAsFactors = FALSE)
    if (verbose) {
        cat("  Binary releases in ", repo, "\n", sep = "")
        for (i in seq_len(nrow(out))) {
            cat(sprintf("  %-22s %-12s %s%s%s\n", out$tag[i],
                        if (is.na(out$published[i])) "" else out$published[i],
                        if (out$installed[i]) "[installed]" else "",
                        if (out$active[i])    " [in use]"   else "",
                        if (out$required[i])  " [required by this package]" else ""))
        }
    }
    invisible(out)
}
