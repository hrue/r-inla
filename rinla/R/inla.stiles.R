#' sTiles helpers: streaming spy plots and graph<->mtx I/O
#'
#' Streaming sparsity-pattern plots and conversions between INLA graph
#' files and MatrixMarket `.mtx` files. The spy plots stream the input
#' file (never loading the full matrix), bin entries into a `GRID x GRID`
#' grid, and render a log-density image. This makes them usable on
#' matrices that would otherwise exhaust memory if read into a `Matrix`.
#'
#' INLA graph files are autodetected as 0-based or 1-based by inspecting
#' the first data row; override with `base = 0L` or `base = 1L`.
#' MatrixMarket files use the standard 1-based coordinate convention.
#'
#' `inla.stiles.graph2mtx` writes a symmetric MatrixMarket file. By
#' default the diagonal is set to the node degree (matches the typical
#' precision-matrix convention where `Q[i,i]` equals the sum of weights of
#' incident edges, giving a diagonally-dominant Q with unit off-diagonals).
#' Use `diagonal = "one"` for unit diagonals, or `diagonal = "none"` to
#' omit the diagonal entirely.
#'
#' @aliases inla.stiles.spy.mtx inla.stiles.spy.graph
#'   inla.stiles.graph2mtx inla.stiles.mtx2graph
#' @param file        For `inla.stiles.spy.mtx`: path to a `.mtx` file.
#'                    For `inla.stiles.spy.graph`: path to an INLA graph file.
#' @param graph_file  Path to an INLA graph file (ASCII).
#' @param mtx_file    Path to a MatrixMarket `.mtx` file.
#' @param out         Output PNG path. Defaults to the input path with
#'                    extension replaced by `.png`.
#' @param GRID        Bin grid resolution (default 1000).
#' @param text        If `TRUE`, draw title and axis labels (default `FALSE`).
#' @param colorbar    If `TRUE`, draw a colorbar to the right (default `FALSE`).
#' @param palette     Colour palette: `"grey"` (default), `"gray"`, `"binary"`,
#'                    `"viridis"`, or any name accepted by `hcl.colors()`.
#' @param diagonal    Diagonal value for `inla.stiles.graph2mtx`:
#'                    `"degree"` (default), `"one"`, or `"none"`.
#' @param base        Index base for INLA graph files. `NULL` (default)
#'                    auto-detects by peeking at the first data row.
#' @param CHUNK       Lines read per streaming chunk. Tune for very large files.
#' @return The output file path (invisibly).
#' @author Esmail Abdul Fattah \email{esmail.abdulfattah@@kaust.edu.sa}
#' @seealso [inla.spy()], [inla.read.graph()], [inla.write.graph()]
#' @examples
#' \dontrun{
#' ## spy a large MTX file (streamed, never fully loaded)
#' inla.stiles.spy.mtx("Q.mtx", "Q.png")
#'
#' ## spy an INLA graph file directly
#' inla.stiles.spy.graph("inla_graph_animal1", "animal.png")
#'
#' ## convert between formats (round-trips bit-identical up to neighbour order)
#' inla.stiles.graph2mtx("inla_graph_animal1", "animal.mtx")
#' inla.stiles.mtx2graph("animal.mtx",         "animal.graph")
#' }
#' @name inla.stiles
NULL

# ---- internal helpers (not exported) ----------------------------------

`inla.stiles.read.mtx.header` <- function(con) {
    hdr <- tolower(trimws(readLines(con, n = 1L, warn = FALSE)))
    is_symmetric <- grepl("symmetric", hdr, fixed = TRUE)
    is_pattern   <- grepl("pattern",   hdr, fixed = TRUE)
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (length(line) == 0L) stop("unexpected EOF before size line")
        if (!startsWith(line, "%")) break
    }
    dims <- as.integer(strsplit(trimws(line), "\\s+")[[1]])
    list(n = dims[1], m = dims[2], nnz = dims[3],
         is_symmetric = is_symmetric, is_pattern = is_pattern)
}

`inla.stiles.detect.graph.base` <- function(file) {
    con <- file(file, open = "r"); on.exit(close(con))
    readLines(con, n = 1L, warn = FALSE)  # skip n
    repeat {
        line <- readLines(con, n = 1L, warn = FALSE)
        if (length(line) == 0L) return(0L)
        line <- trimws(line)
        if (nchar(line) > 0L) {
            first <- as.integer(strsplit(line, "\\s+")[[1]][1])
            return(if (!is.na(first) && first == 0L) 0L else 1L)
        }
    }
}

`inla.stiles.add.to.bins` <- function(grid, GRID, ii, jj, mirror) {
    idx <- ii + jj * GRID + 1L
    grid <- grid + matrix(tabulate(idx, nbins = GRID * GRID),
                          nrow = GRID, ncol = GRID)
    if (length(mirror) > 0L) {
        idx2 <- jj[mirror] + ii[mirror] * GRID + 1L
        grid <- grid + matrix(tabulate(idx2, nbins = GRID * GRID),
                              nrow = GRID, ncol = GRID)
    }
    grid
}

`inla.stiles.render.spy` <- function(grid, n, nnz, title, out_path, GRID,
                                     text = FALSE, colorbar = FALSE,
                                     palette = "grey") {
    img <- ifelse(grid > 0, log10(grid + 1), NA_real_)
    pal <- switch(palette,
        grey    = colorRampPalette(c("#f0f0f0", "black"))(256),
        gray    = colorRampPalette(c("#f0f0f0", "black"))(256),
        binary  = colorRampPalette(c("white",   "black"))(256),
        viridis = hcl.colors(256, palette = "viridis"),
        hcl.colors(256, palette = palette))
    zmax <- suppressWarnings(max(img, na.rm = TRUE))
    if (!is.finite(zmax)) zmax <- 1
    scale <- GRID / n

    png(out_path,
        width  = if (colorbar) 2700 else 2400,
        height = 2400, res = 300,
        type = if (capabilities("cairo")) "cairo-png" else "Xlib")
    op <- par(no.readonly = TRUE); on.exit(par(op), add = TRUE)
    if (colorbar) layout(matrix(c(1, 2), nrow = 1), widths = c(6, 1))

    par(mar = if (text) c(4, 4.5, 3.5, 1) else c(1, 1, 1, 1))
    image(x = seq_len(GRID), y = seq_len(GRID),
          z = t(img)[, GRID:1, drop = FALSE],
          col = pal, zlim = c(0, zmax),
          xlab = if (text) "column" else "",
          ylab = if (text) "row"    else "",
          axes = FALSE, useRaster = TRUE,
          main = if (text)
                     sprintf("%s\nn=%s  nnz=%s  grid=%dx%d",
                             title, format(n, big.mark = ","),
                             format(nnz, big.mark = ","), GRID, GRID)
                 else "")
    if (text) {
        ax <- pretty(c(1, GRID))
        axis(1, at = ax, labels = format(round(ax / scale), big.mark = ","))
        axis(2, at = ax,
             labels = format(round((GRID - ax) / scale), big.mark = ","),
             las = 1)
    }
    box()

    if (colorbar) {
        par(mar = if (text) c(4, 0.5, 3.5, 4) else c(3, 0.5, 1, 3))
        zs <- seq(0, zmax, length.out = 256)
        image(x = 1, y = zs, z = matrix(zs, nrow = 1),
              col = pal, axes = FALSE, xlab = "", ylab = "", useRaster = TRUE)
        axis(4, las = 1)
        if (text) mtext("log10(nnz per cell + 1)", side = 4, line = 2.6,
                        cex = 0.85)
        box()
    }
    dev.off()
    message(sprintf("wrote %s", out_path))
}

# ---- public API -------------------------------------------------------

#' @rdname inla.stiles
#' @export
`inla.stiles.spy.mtx` <- function(file, out = NULL, GRID = 1000L,
                                  text = FALSE, colorbar = FALSE,
                                  palette = "grey",
                                  CHUNK = 2000000L) {
    if (is.null(out)) out <- sub("\\.[^.]*$", ".png", file)
    con <- file(file, open = "r"); on.exit(close(con), add = TRUE)
    h <- inla.stiles.read.mtx.header(con)
    stopifnot(h$n == h$m)

    scale <- GRID / h$n
    grid  <- matrix(0L, GRID, GRID)
    count <- 0L
    repeat {
        lines <- readLines(con, n = CHUNK, warn = FALSE)
        if (length(lines) == 0L) break
        parts <- do.call(rbind, strsplit(lines, "\\s+"))
        i <- as.integer(parts[, 1]) - 1L
        j <- as.integer(parts[, 2]) - 1L
        keep <- !is.na(i) & !is.na(j)
        i <- i[keep]; j <- j[keep]
        ii <- pmin(as.integer(i * scale), GRID - 1L)
        jj <- pmin(as.integer(j * scale), GRID - 1L)
        mirror <- if (h$is_symmetric) which(i != j) else integer(0)
        grid  <- inla.stiles.add.to.bins(grid, GRID, ii, jj, mirror)
        count <- count + length(i)
    }
    message(sprintf("binned %s entries into %dx%d",
                    format(count, big.mark = ","), GRID, GRID))
    inla.stiles.render.spy(grid, h$n, h$nnz, basename(file), out, GRID,
                           text = text, colorbar = colorbar, palette = palette)
    invisible(out)
}

#' @rdname inla.stiles
#' @export
`inla.stiles.spy.graph` <- function(file, out = NULL, GRID = 1000L,
                                    text = FALSE, colorbar = FALSE,
                                    palette = "grey",
                                    CHUNK = 100000L, base = NULL) {
    if (is.null(out))  out  <- sub("\\.[^.]*$", ".png", file)
    if (is.null(base)) base <- inla.stiles.detect.graph.base(file)

    con <- file(file, open = "r"); on.exit(close(con), add = TRUE)
    n <- as.integer(trimws(readLines(con, n = 1L, warn = FALSE)))
    stopifnot(!is.na(n), n > 0L)

    scale <- GRID / n
    grid  <- matrix(0L, GRID, GRID)
    nnz_off <- 0L
    repeat {
        lines <- readLines(con, n = CHUNK, warn = FALSE)
        if (length(lines) == 0L) break
        toks <- strsplit(trimws(lines), "\\s+")
        toks <- toks[lengths(toks) >= 2L]
        if (length(toks) == 0L) next

        nb <- vapply(toks, function(t) as.integer(t[2]), integer(1))
        rep_i <- rep.int(vapply(toks, function(t) as.integer(t[1]),
                                integer(1)), nb)
        nbrs  <- unlist(lapply(toks, function(t) {
            k <- as.integer(t[2])
            if (k > 0L) as.integer(t[3:(2 + k)]) else integer(0)
        }), use.names = FALSE)

        i0 <- rep_i - base
        j0 <- nbrs  - base
        ii <- pmin(as.integer(i0 * scale), GRID - 1L)
        jj <- pmin(as.integer(j0 * scale), GRID - 1L)
        grid <- inla.stiles.add.to.bins(grid, GRID, ii, jj, integer(0))
        nnz_off <- nnz_off + length(i0)
    }
    diag_idx <- pmin(as.integer(seq.int(0L, n - 1L) * scale), GRID - 1L)
    for (k in diag_idx) grid[k + 1L, k + 1L] <- grid[k + 1L, k + 1L] + 1L

    nnz_total <- nnz_off + n
    message(sprintf("graph: n=%s  base=%d  nnz(incl. diag, both dirs)=%s",
                    format(n, big.mark = ","), base,
                    format(nnz_total, big.mark = ",")))
    inla.stiles.render.spy(grid, n, nnz_total, basename(file), out, GRID,
                           text = text, colorbar = colorbar, palette = palette)
    invisible(out)
}

#' @rdname inla.stiles
#' @export
`inla.stiles.graph2mtx` <- function(graph_file, mtx_file,
                                    diagonal = c("degree", "one", "none"),
                                    base = NULL) {
    diagonal <- match.arg(diagonal)
    if (is.null(base)) base <- inla.stiles.detect.graph.base(graph_file)

    lines <- readLines(graph_file, warn = FALSE)
    n <- as.integer(trimws(lines[1])); stopifnot(!is.na(n), n > 0L)
    toks <- strsplit(trimws(lines[-1]), "\\s+")
    toks <- toks[lengths(toks) >= 2L]

    row_idx <- vapply(toks, function(t) as.integer(t[1]), integer(1))
    degs    <- vapply(toks, function(t) as.integer(t[2]), integer(1))
    nbrs    <- lapply(toks, function(t) {
        k <- as.integer(t[2])
        if (k > 0L) as.integer(t[3:(2 + k)]) else integer(0)
    })

    i_all <- rep.int(row_idx - base + 1L, degs)
    j_all <- unlist(nbrs, use.names = FALSE) - base + 1L
    upper <- i_all < j_all
    i_up <- i_all[upper]; j_up <- j_all[upper]

    deg_by_row <- integer(n)
    deg_by_row[row_idx - base + 1L] <- degs

    out_con <- file(mtx_file, open = "w")
    on.exit(close(out_con), add = TRUE)
    writeLines("%%MatrixMarket matrix coordinate real symmetric", out_con)

    if (diagonal == "none") {
        writeLines(sprintf("%d %d %d", n, n, length(i_up)), out_con)
    } else {
        diag_v <- if (diagonal == "degree") deg_by_row else rep_len(1.0, n)
        writeLines(sprintf("%d %d %d", n, n, length(i_up) + n), out_con)
        writeLines(sprintf("%d %d %g", seq_len(n), seq_len(n),
                           as.numeric(diag_v)), out_con)
    }
    if (length(i_up) > 0L)
        writeLines(sprintf("%d %d 1.0", i_up, j_up), out_con)

    message(sprintf("wrote %s  (n=%d, nnz_off_upper=%d, diagonal=%s, base=%d)",
                    mtx_file, n, length(i_up), diagonal, base))
    invisible(mtx_file)
}

#' @rdname inla.stiles
#' @export
`inla.stiles.mtx2graph` <- function(mtx_file, graph_file,
                                    base = 0L, CHUNK = 2000000L) {
    con <- file(mtx_file, open = "r"); on.exit(close(con), add = TRUE)
    h <- inla.stiles.read.mtx.header(con)
    stopifnot(h$n == h$m)
    n <- h$n

    acc_i <- vector("list", 0L); acc_j <- vector("list", 0L); k <- 0L
    repeat {
        lines <- readLines(con, n = CHUNK, warn = FALSE)
        if (length(lines) == 0L) break
        parts <- do.call(rbind, strsplit(lines, "\\s+"))
        i <- as.integer(parts[, 1]); j <- as.integer(parts[, 2])
        keep <- !is.na(i) & !is.na(j) & i != j
        i <- i[keep]; j <- j[keep]
        if (h$is_symmetric) { a <- c(i, j); b <- c(j, i) } else { a <- i; b <- j }
        k <- k + 1L; acc_i[[k]] <- a; acc_j[[k]] <- b
    }
    a <- unlist(acc_i, use.names = FALSE)
    b <- unlist(acc_j, use.names = FALSE)

    key <- (as.numeric(a) - 1) * n + as.numeric(b)
    dup <- duplicated(key)
    a <- a[!dup]; b <- b[!dup]

    ord <- order(a, b)
    a <- a[ord]; b <- b[ord]
    nbrs_by_row <- split(b, factor(a, levels = seq_len(n)))

    out <- file(graph_file, open = "w"); on.exit(close(out), add = TRUE)
    writeLines(as.character(n), out)
    buf <- character(n)
    for (i in seq_len(n)) {
        nb <- nbrs_by_row[[i]]
        nb_out <- nb - 1L + base
        i_out  <- i  - 1L + base
        if (length(nb) == 0L) {
            buf[i] <- sprintf("%d 0", i_out)
        } else {
            buf[i] <- paste(i_out, length(nb), paste(nb_out, collapse = " "))
        }
    }
    writeLines(buf, out)
    message(sprintf("wrote %s  (n=%d, base=%d)", graph_file, n, base))
    invisible(graph_file)
}
