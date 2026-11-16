# Shared plotting helpers: theme, palette, devices, genome layout, annotation.

#' House ggplot theme
#'
#' A modification of [ggplot2::theme_linedraw()]: no grid lines, transparent strip background with black
#' strip text, legend at the bottom, base font size 6. Used by every `rae_plot_*` function.
#'
#' @param base_size base font size in points (default 6).
#' @param base_family font family; `""` is the device default.
#' @param ... further arguments passed to [ggplot2::theme()].
#' @return a ggplot2 theme.
#' @export
theme_rodrigo <- function(base_size = 6, base_family = "", ...) {
  .need("ggplot2")
  ggplot2::`%+replace%`(
    ggplot2::theme_linedraw(base_size = base_size, base_family = base_family),
    ggplot2::theme(panel.grid = ggplot2::element_blank(),
                   legend.position = "bottom",
                   strip.background = ggplot2::element_rect(colour = "#FFFFFF00", fill = "#FFFFFF00"),
                   strip.text = ggplot2::element_text(colour = "#000000"),
                   ...))
}

# Palette. Gain/loss pass the colour-vision-deficiency checks against each other and stay distinct in
# greyscale; germline-overlapping events are grey and carry the mark "Г".
rae_colours <- c(gain = "#D40000", loss = "#4A90D9", germline = "#8F8F8F", ns = "#BABABA",
                 ink = "#000000", ink2 = "#404040")

.need <- function(pkgs) {
  miss <- pkgs[!vapply(pkgs, requireNamespace, logical(1), quietly = TRUE)]
  if (length(miss)) stop("package(s) needed for this plot: ", paste(miss, collapse = ", "), call. = FALSE)
}

#' Save a figure with the base graphics devices
#'
#' Opens [grDevices::png()] or [grDevices::cairo_pdf()] (chosen from the file extension) at the given size
#' and draws `p` (a ggplot or a ComplexHeatmap object). `cairo_pdf()` is used rather than `pdf()` because
#' `pdf()` with the default fonts cannot draw the germline mark "Г".
#'
#' @param p a ggplot, a ComplexHeatmap `Heatmap`/`HeatmapList`, or a function drawing on the device.
#' @param file output path ending in `.png` or `.pdf`.
#' @param width,height size in millimetres.
#' @param res resolution for PNG, pixels per inch.
#' @param ... passed to `ComplexHeatmap::draw()` for heatmaps.
#' @return `file`, invisibly.
#' @export
rae_save <- function(p, file, width = 100, height = 100, res = 600, ...) {
  ext <- tolower(tools::file_ext(file))
  # PNG device: on macOS the quartz device. A cairo PNG at 600 dpi and cairo_pdf() in the same R session share
  # font metrics across resolutions, and every later figure is laid out with wrong text sizes (titles over tick
  # labels, swollen legends); quartz does not share them. Quartz writes 72 dpi into the file, so the resolution
  # tag is rewritten afterwards (package 'png', if installed). Both devices draw "Г".
  quartz <- ext == "png" && isTRUE(capabilities("aqua"))
  if (ext == "png") {
    grDevices::png(file, width = width, height = height, units = "mm", res = res,
                   type = if (quartz) "quartz" else if (capabilities("cairo")) "cairo" else getOption("bitmapType"))
    if (quartz && requireNamespace("png", quietly = TRUE))
      on.exit(png::writePNG(png::readPNG(file), file, dpi = res), add = TRUE, after = TRUE)
  } else if (ext == "pdf") grDevices::cairo_pdf(file, width = width / 25.4, height = height / 25.4)
  else stop("file must end in .png or .pdf")
  on.exit(grDevices::dev.off(), add = TRUE, after = FALSE)
  if (inherits(p, c("Heatmap", "HeatmapList"))) ComplexHeatmap::draw(p, ...)
  else if (is.function(p)) p()
  else print(p)
  invisible(file)
}

#' Cytoband table of a genome build
#'
#' @param build `"hg38"`, `"hg19"` or `"chm13"` (aliases `"hs1"`, `"t2t"`). Source: UCSC (free for any use).
#' @return data frame with chrom (`"1"`..`"22"`, `"X"`, `"Y"`), start (1-based), end, band, stain; `NULL` for
#'   other builds.
#' @export
rae_cytobands <- function(build = "hg38") {
  b <- c(hs1 = "chm13", t2t = "chm13", chm13 = "chm13", hg19 = "hg19", hg38 = "hg38")[tolower(build)]
  if (is.na(b)) return(NULL)
  f <- system.file("extdata", paste0("cytoband_", b, ".tsv.gz"), package = "rae3")
  if (!nzchar(f)) return(NULL)
  utils::read.delim(f, colClasses = c("character", "numeric", "numeric", "character", "character"))
}

# ---- genome layout ----------------------------------------------------------------------------------

# Chromosome lengths and centromeres for the result's build (cytobands, else the core's annotation, else data).
.chrom_info <- function(x) {
  ch <- sort(unique(as.integer(x$regions$chrom)))
  dmax <- tapply(x$regions$end, x$regions$chrom, max)
  d <- x$diagnostics
  len <- if (!is.null(d$chrom_len) && length(d$chrom_len) >= max(ch) && all(d$chrom_len[ch] > 0)) d$chrom_len[ch] else as.numeric(dmax[as.character(ch)])
  cb <- rae_cytobands(if (is.null(x$build)) "custom" else x$build)
  cs <- ce <- rep(NA_real_, length(ch))
  if (!is.null(cb)) {
    a <- cb[cb$stain == "acen", ]
    cs <- as.numeric(tapply(a$start, a$chrom, min)[as.character(ch)])
    ce <- as.numeric(tapply(a$end, a$chrom, max)[as.character(ch)])
  } else if (!is.null(d$cen_start) && length(d$cen_start) >= max(ch)) {
    cs <- ifelse(d$cen_start[ch] > 0, d$cen_start[ch], NA); ce <- ifelse(d$cen_end[ch] > 0, d$cen_end[ch], NA)
  }
  data.frame(chrom = ch, len = pmax(len, as.numeric(dmax[as.character(ch)])), cen_start = cs, cen_end = ce)
}

# Units to draw, in order: whole chromosomes or arms, from a spec such as c(8, "10q", "chr13").
.layout <- function(x, chrom = NULL, gap = 0) {
  ci <- .chrom_info(x)
  spec <- if (is.null(chrom)) as.character(ci$chrom) else sub("^chr", "", as.character(chrom), ignore.case = TRUE)
  u <- lapply(spec, function(s) {
    arm <- if (grepl("[pq]$", s)) substring(s, nchar(s)) else ""
    c0 <- as.integer(sub("[pq]$", "", s))
    i <- match(c0, ci$chrom)
    if (is.na(i)) stop("chromosome not in the result: ", s, call. = FALSE)
    if (nzchar(arm) && is.na(ci$cen_start[i])) stop("centromere unknown for chromosome ", c0, "; arms cannot be selected", call. = FALSE)
    from <- if (arm == "q") ci$cen_end[i] + 1 else 1
    to <- if (arm == "p") ci$cen_start[i] - 1 else ci$len[i]
    data.frame(chrom = c0, arm = arm, from = from, to = to, cen_start = ci$cen_start[i], cen_end = ci$cen_end[i],
               label = paste0(c0, arm))
  })
  u <- do.call(rbind, u)
  w <- u$to - u$from + 1
  u$offset <- c(0, cumsum(w + gap))[seq_len(nrow(u))]
  u
}

# Genome coordinate of (chrom, pos) in the layout; NA outside every unit. Intervals are clipped to units.
.gpos <- function(L, chrom, pos) {
  out <- rep(NA_real_, length(pos))
  for (k in seq_len(nrow(L))) {
    i <- chrom == L$chrom[k] & pos >= L$from[k] & pos <= L$to[k]
    out[i] <- L$offset[k] + pos[i] - L$from[k]
  }
  out
}
.ginterval <- function(L, chrom, start, end) {
  res <- vector("list", nrow(L))
  for (k in seq_len(nrow(L))) {
    i <- which(chrom == L$chrom[k] & end >= L$from[k] & start <= L$to[k])
    if (!length(i)) next
    res[[k]] <- data.frame(row = i, xmin = L$offset[k] + pmax(start[i], L$from[k]) - L$from[k],
                           xmax = L$offset[k] + pmin(end[i], L$to[k]) - L$from[k])
  }
  do.call(rbind, res)
}

# ---- labels: cytobands and cancer genes ---------------------------------------------------------------

.band_label <- function(cb, chrom, start, end) {
  if (is.null(cb)) return(paste0(chrom, ":", format(round(start / 1e6, 2), nsmall = 2), "-", format(round(end / 1e6, 2), nsmall = 2), " Mb"))
  b <- cb[cb$chrom == as.character(chrom) & cb$end >= start & cb$start <= end, "band"]
  if (!length(b)) return(as.character(chrom))
  if (length(b) == 1) b else paste0(b[1], "–", sub("^[0-9XY]+", "", b[length(b)]))
}

.genes_table <- function(genes, build) {
  if (is.null(genes)) return(NULL)
  if (is.data.frame(genes)) return(genes)
  if (identical(genes, "oncokb")) {
    g <- tryCatch(rae_cancer_genes(build), error = function(e) NULL)
    if (is.null(g)) message("no cancer-gene list for build '", build, "': labelling by cytoband only. Run rae_fetch_cancer_genes(\"",
                            build, "\") once to add gene names.")
    return(g)
  }
  stop("genes must be \"oncokb\", a data frame (gene, chrom, start, end) or NULL")
}

# One label per event: germline mark, focal cancer genes, or cytoband(s).
.event_label <- function(ev, cb, gt, focal_max_bp, focal_max_genes, symbol) {
  vapply(seq_len(nrow(ev)), function(i) {
    if (isTRUE(ev$germline_overlap[i])) return(symbol)
    if (!is.null(gt) && ev$size_bp[i] <= focal_max_bp) {
      g <- gt[as.character(gt$chrom) == as.character(ev$chrom[i]) & gt$end >= ev$start[i] & gt$start <= ev$end[i], "gene"]
      if (length(g) >= 1 && length(g) <= focal_max_genes) return(paste(unique(g), collapse = "\n"))
    }
    .band_label(cb, ev$chrom[i], ev$start[i], ev$end[i])
  }, character(1))
}

# Events to label: significant peaks, plus significant stage-1 events that hold no peak.
.label_events <- function(x, fdr) {
  r <- x$roi
  if (is.null(r$germline_overlap)) r$germline_overlap <- FALSE
  r <- r[r$q <= fdr, , drop = FALSE]
  pk <- r[r$type == "peak", , drop = FALSE]
  ev <- r[r$type == "region", , drop = FALSE]
  has_peak <- vapply(seq_len(nrow(ev)), function(i)
    any(pk$direction == ev$direction[i] & pk$chrom == ev$chrom[i] & pk$start >= ev$start[i] & pk$end <= ev$end[i]), logical(1))
  out <- rbind(pk, ev[!has_peak, , drop = FALSE])
  out[order(out$q, -out$score), , drop = FALSE]
}

# Germline intervals per direction (kept and flagged, or filtered out).
.germline_intervals <- function(x) {
  g <- rbind(if (!is.null(x$roi$germline_overlap)) x$roi[x$roi$germline_overlap, , drop = FALSE], x$roi_germline)
  if (is.null(g)) g <- x$roi[0, , drop = FALSE]
  g[, c("direction", "chrom", "start", "end")]
}
.in_intervals <- function(chrom, start, end, iv) {
  hit <- rep(FALSE, length(start))
  for (k in seq_len(nrow(iv))) hit <- hit | (chrom == iv$chrom[k] & start <= iv$end[k] & end >= iv$start[k])
  hit
}

# Score at which q first reaches each level (q is non-increasing in the score), for a q-labelled axis.
.q_axis <- function(score, q, fdr, gap = NULL) {
  ok <- is.finite(score) & is.finite(q)
  score <- score[ok]; q <- q[ok]
  if (!length(q)) return(list(at = numeric(), lab = character()))
  lev <- c(fdr, 10^-(2:12))
  lev <- lev[lev >= min(q)]
  at <- vapply(lev, function(l) min(score[q <= l]), numeric(1))
  lab <- ifelse(lev == fdr, format(fdr), paste0("10^", log10(lev)))
  if (any(q <= min(q) * (1 + 1e-9)) && min(q) < fdr) {
    at <- c(at, max(score)); lab <- c(lab, paste0("'≤'*", "10^", round(log10(min(q)), 1)))
  }
  # thin ticks that would collide: keep the first (the FDR level) and the last (the floor)
  keep <- rep(FALSE, length(at)); last <- -Inf; if (is.null(gap)) gap <- 0.06 * max(score)
  for (i in seq_along(at)) if (i == 1 || at[i] - last >= gap) { keep[i] <- TRUE; last <- at[i] }
  n <- length(at)
  if (n > 1 && !keep[n]) { keep[max(which(keep))] <- max(which(keep)) == 1; keep[n] <- TRUE }
  list(at = at[keep], lab = parse(text = lab[keep]))
}
