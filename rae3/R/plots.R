# Default figures. ggplot2 figures return a ggplot; the heatmap returns a ComplexHeatmap object.
# Save with rae_save() (png() / cairo_pdf(); never ggsave). Copy number is log2 ratio throughout.

.dir_col <- function(d) unname(rae_colours[ifelse(d == "gain", "gain", "loss")])

#' Genome-wide significant alterations (Taylor 2010 Fig. 1A; Taylor 2008 Fig. 5)
#'
#' Summary scores of the unified breakpoint profile, A' (gain/amplification) upwards and D' (loss/deletion)
#' downwards, with the chromosome ideogram in the centre. Right axis: score; left axis: the FDR (q) reached at
#' that score. Each significant peak (and each significant event without a peak) is labelled with its cancer
#' genes when it is focal (at most `focal_max_bp` long and holding 1 to `focal_max_genes` listed genes),
#' otherwise with its cytoband(s). Events overlapping germline variation are drawn grey and marked with
#' `germline_symbol`; they are absent if the result was made with `filter_germline_overlaps = TRUE`.
#'
#' @param x a `rae_result`.
#' @param chrom `NULL` for the genome, or chromosomes and/or arms to draw, e.g. `c(8, "10q", "13q")`.
#' @param genes `"oncokb"` (the locally fetched list, see [rae_fetch_cancer_genes()]), a data frame with
#'   columns gene, chrom, start, end, or `NULL` for cytoband labels only.
#' @param focal_max_bp,focal_max_genes focal-label rule (defaults 3 Mb and 10 genes).
#' @param show `"significant"` (Taylor 2010: only regions with q <= `fdr`) or `"all"` (Taylor 2008: all regions,
#'   the non-significant ones in light grey).
#' @param fdr significance threshold; defaults to the run's `fdr`.
#' @param max_labels most labels per direction, best q first.
#' @param germline_symbol mark for germline overlaps.
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_genome <- function(x, chrom = NULL, genes = "oncokb", focal_max_bp = 3e6, focal_max_genes = 10,
                            show = c("significant", "all"), fdr = x$params$fdr, max_labels = 15,
                            germline_symbol = "Г", base_size = 6) {
  .need(c("ggplot2", "ggrepel"))
  show <- match.arg(show)
  L <- .layout(x, chrom)
  r <- x$regions
  smax <- max(c(r$A_score, r$D_score), na.rm = TRUE)
  h <- 0.045 * smax                                  # half-height of the ideogram
  giv <- .germline_intervals(x)
  bars <- list()
  for (d in c("gain", "loss")) {
    sc <- if (d == "gain") r$A_score else r$D_score
    q <- if (d == "gain") r$A_q else r$D_q
    iv <- .ginterval(L, r$chrom, r$start, r$end)
    if (is.null(iv)) next
    sig <- q[iv$row] <= fdr
    germ <- .in_intervals(r$chrom[iv$row], r$start[iv$row], r$end[iv$row], giv[giv$direction == d, ])
    col <- ifelse(germ, "germline", ifelse(sig, d, "ns"))
    keep <- if (show == "significant") sig else rep(TRUE, length(sig))
    if (isTRUE(x$filter_germline_overlaps)) keep <- keep & !germ
    s <- sc[iv$row]
    bars[[d]] <- data.frame(xmin = iv$xmin, xmax = iv$xmax, ymin = if (d == "gain") h else -h,
                            ymax = if (d == "gain") h + s else -h - s, col = col)[keep, , drop = FALSE]
  }
  bars <- do.call(rbind, bars)
  # ideogram
  # chromosomes alternate black and white; the score-0 lines along both edges close the white ones
  odd <- seq_len(nrow(L)) %% 2 == 1
  ide <- data.frame(xmin = L$offset, xmax = L$offset + L$to - L$from, ymin = -h, ymax = h,
                    fill = ifelse(odd, "#000000", "#FFFFFF"))
  cen <- data.frame(x = .gpos(L, L$chrom, (L$cen_start + L$cen_end) / 2), col = ifelse(odd, "#FFFFFF", "#000000"))
  cen <- cen[is.finite(cen$x), , drop = FALSE]
  lab <- data.frame(x = (ide$xmin + ide$xmax) / 2, y = 0, label = L$label, col = ifelse(odd, "#FFFFFF", "#000000"))
  # axes: score on the right, q on the left
  qa <- .q_axis(r$A_score, r$A_q, fdr, gap = 0.1 * smax); qd <- .q_axis(r$D_score, r$D_q, fdr, gap = 0.1 * smax)
  br <- pretty(c(0, smax), 4); br <- br[br > 0 & br <= smax * 1.05]
  thr <- c(if (length(qa$at)) h + qa$at[1], if (length(qd$at)) -h - qd$at[1])
  # the same score range above and below the ideogram, with room for two label rows beyond the tallest bar
  ymax <- function(v) if (length(v)) max(v, na.rm = TRUE) else 0
  tallest <- max(ymax(r$D_score[r$D_q <= fdr | show == "all"]), ymax(r$A_score[r$A_q <= fdr | show == "all"]), 0.1 * smax)
  step <- 0.12 * tallest
  m <- h + tallest + 0.03 * smax + 2 * step + 0.1 * tallest
  lim <- c(-m, m)
  # labels
  ev <- .label_events(x, fdr)
  ev <- ev[!is.na(.gpos(L, ev$chrom, (ev$start + ev$end) / 2)), , drop = FALSE]
  ev <- do.call(rbind, lapply(split(ev, ev$direction), function(e) utils::head(e, max_labels)))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = ide, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax),
                       fill = ide$fill, colour = NA)
  if (nrow(cen)) {
    cs <- rbind(transform(cen, y0 = 0.6 * h, y1 = h), transform(cen, y0 = -h, y1 = -0.6 * h))
    p <- p + ggplot2::geom_segment(data = cs, ggplot2::aes(x = x, xend = x, y = y0, yend = y1), colour = cs$col, linewidth = 0.4)
  }
  p <- p + ggplot2::geom_hline(yintercept = c(-h, h), linewidth = 0.25, colour = rae_colours[["ink"]])
  if (!is.null(bars) && nrow(bars))
    p <- p + ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax,
                                                          fill = col, colour = col), linewidth = 0.15) +
      ggplot2::scale_fill_manual(values = rae_colours, breaks = c("gain", "loss", "germline"),
                                 labels = c("gain / amplification (A')", "loss / deletion (D')",
                                            paste("overlaps germline variation", germline_symbol)), name = NULL,
                                 aesthetics = c("fill", "colour"))
  p <- p + ggplot2::geom_text(data = lab, ggplot2::aes(x = x, y = y, label = label), colour = lab$col,
                              size = base_size * 0.7 / ggplot2::.pt) +
    ggplot2::geom_hline(yintercept = thr, linetype = "22", linewidth = 0.25, colour = rae_colours[["ink2"]])
  if (!is.null(ev) && nrow(ev)) {
    gt <- .genes_table(genes, x$build)
    ev$label <- .event_label(ev, rae_cytobands(x$build), gt, focal_max_bp, focal_max_genes, germline_symbol)
    ev$x <- .gpos(L, ev$chrom, (ev$start + ev$end) / 2)
    ev$y <- ifelse(ev$direction == "gain", h + ev$score, -h - ev$score)
    ev$face <- ifelse(ev$germline_overlap | grepl("^[0-9XY]+[pq]|:", ev$label), "plain", "italic")
    for (d in c("gain", "loss")) {
      e <- ev[ev$direction == d, , drop = FALSE]
      if (!nrow(e)) next
      sgn <- if (d == "gain") 1 else -1
      bd <- bars[bars$ymin == sgn * h, , drop = FALSE]
      # two label rows beyond the tallest bar of this direction; labels alternate rows in genome order and
      # are spread sideways only, so they never sit in the data
      top <- if (nrow(bd)) max(abs(bd$ymax)) else h
      e <- e[order(e$x), , drop = FALSE]
      e$row <- (seq_len(nrow(e)) - 1) %% 2
      e$nud <- sgn * (top + 0.03 * smax + e$row * step) - e$y
      e$col <- ifelse(e$germline_overlap, "germline", "ink2")
      p <- p + ggrepel::geom_text_repel(
        data = e, ggplot2::aes(x = x, y = y, label = label, fontface = face), size = base_size / ggplot2::.pt,
        colour = unname(rae_colours[e$col]), lineheight = 0.85, direction = "x", segment.size = 0.15,
        segment.colour = rae_colours[["ns"]], min.segment.length = 0, box.padding = 0.12, max.overlaps = Inf,
        nudge_y = e$nud, vjust = if (d == "gain") 0 else 1, seed = 1)
    }
  }
  p + ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(
      name = "score (A' above, D' below)", position = "right", limits = lim,
      breaks = c(h + br, -h - br), labels = format(c(br, br)),
      sec.axis = ggplot2::sec_axis(~ ., name = "q (FDR)", breaks = c(h + qa$at, -h - qd$at),
                                   labels = c(qa$lab, qd$lab))) +
    theme_rodrigo(base_size) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_line(linewidth = 0.25),
                   legend.key.size = ggplot2::unit(3, "mm"))
}

#' Genome-wide alteration frequency
#'
#' Percentage of tumours called in each state, per region of the unified breakpoint profile (the discrete
#' calls of `calls.tsv`: `call_threshold`, `amp_threshold`), gains and amplifications upwards and losses and
#' deletions downwards, with the chromosome ideogram in the centre. Unlike [rae_plot_genome()], which shows the
#' summary score and the FDR it reaches, this shows only how often each event occurs; it says nothing about
#' whether that is more than expected by chance.
#'
#' @param x a `rae_result`.
#' @param type `"stacked"` (default): single-copy gain with amplification stacked on top of it above the
#'   ideogram, single-copy loss with deletion stacked beneath it below. Otherwise one or more of `"gain"`
#'   (all gains, amplification included), `"amplification"`, `"loss"` (all losses, deletion included) and
#'   `"deletion"`; gains and amplifications are drawn upwards, losses and deletions downwards.
#' @param chrom `NULL` for the genome, or chromosomes and/or arms, as in [rae_plot_genome()].
#' @param max_pct top of both frequency axes (percent, default 100).
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_frequency <- function(x, type = "stacked", chrom = NULL, max_pct = 100, base_size = 6) {
  .need("ggplot2")
  ok <- c("stacked", "gain", "amplification", "loss", "deletion")
  if (!length(type) || !all(type %in% ok)) stop("type must be \"stacked\" or any of ", paste(dQuote(ok[-1], FALSE), collapse = ", "))
  if ("stacked" %in% type && length(type) > 1) stop("\"stacked\" cannot be combined with other types")
  L <- .layout(x, chrom)
  r <- x$regions
  iv <- .ginterval(L, r$chrom, r$start, r$end)
  if (is.null(iv)) stop("no regions in the chosen chromosomes")
  f <- r[iv$row, c("freq_gain", "freq_amp", "freq_loss", "freq_homdel")] * 100
  h <- 0.045 * max_pct                               # half-height of the ideogram, as in rae_plot_genome()
  # each layer: frequency range [lo, hi] (percent) on one side of the ideogram
  layer <- function(lo, hi, side, state) data.frame(xmin = iv$xmin, xmax = iv$xmax, ymin = side * (h + lo),
                                                    ymax = side * (h + hi), state = state)
  bars <- if (identical(type, "stacked")) {
    rbind(layer(0, f$freq_gain - f$freq_amp, 1, "gain"), layer(f$freq_gain - f$freq_amp, f$freq_gain, 1, "amplification"),
          layer(0, f$freq_loss - f$freq_homdel, -1, "loss"), layer(f$freq_loss - f$freq_homdel, f$freq_loss, -1, "deletion"))
  } else {
    col <- c(gain = "freq_gain", amplification = "freq_amp", loss = "freq_loss", deletion = "freq_homdel")
    if (all(c("gain", "amplification") %in% type) || all(c("loss", "deletion") %in% type))
      stop("gain and amplification (or loss and deletion) overlap; use type = \"stacked\" to show both")
    do.call(rbind, lapply(type, function(t) layer(0, f[[col[[t]]]], if (t %in% c("gain", "amplification")) 1 else -1, t)))
  }
  bars <- bars[is.finite(bars$ymax) & abs(bars$ymax) > abs(bars$ymin), , drop = FALSE]
  bars$state <- factor(bars$state, c("amplification", "gain", "loss", "deletion"))
  pal <- c(amplification = "#8C0000", gain = rae_colours[["gain"]], loss = rae_colours[["loss"]], deletion = "#1B4F8A")
  odd <- seq_len(nrow(L)) %% 2 == 1
  ide <- data.frame(xmin = L$offset, xmax = L$offset + L$to - L$from, ymin = -h, ymax = h,
                    fill = ifelse(odd, "#000000", "#FFFFFF"))
  cen <- data.frame(x = .gpos(L, L$chrom, (L$cen_start + L$cen_end) / 2), col = ifelse(odd, "#FFFFFF", "#000000"))
  cen <- cen[is.finite(cen$x), , drop = FALSE]
  lab <- data.frame(x = (ide$xmin + ide$xmax) / 2, y = 0, label = L$label, col = ifelse(odd, "#FFFFFF", "#000000"))
  up <- any(type %in% c("stacked", "gain", "amplification")); down <- any(type %in% c("stacked", "loss", "deletion"))
  br <- pretty(c(0, max_pct), 4); br <- br[br > 0 & br <= max_pct]
  lim <- c(if (down) -(h + max_pct) else -h, if (up) h + max_pct else h)
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = ide, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax), fill = ide$fill, colour = NA)
  if (nrow(cen)) {
    cs <- rbind(transform(cen, y0 = 0.6 * h, y1 = h), transform(cen, y0 = -h, y1 = -0.6 * h))
    p <- p + ggplot2::geom_segment(data = cs, ggplot2::aes(x = x, xend = x, y = y0, yend = y1), colour = cs$col, linewidth = 0.4)
  }
  p + ggplot2::geom_hline(yintercept = c(-h, h), linewidth = 0.25, colour = rae_colours[["ink"]]) +
    ggplot2::geom_rect(data = bars, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = ymin, ymax = ymax, fill = state,
                                                 colour = state), linewidth = 0.15) +
    ggplot2::scale_fill_manual(values = pal, name = NULL, drop = TRUE, aesthetics = c("fill", "colour")) +
    ggplot2::geom_text(data = lab, ggplot2::aes(x = x, y = y, label = label), colour = lab$col,
                       size = base_size * 0.7 / ggplot2::.pt) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::scale_y_continuous(name = "tumours (%)", limits = lim, expand = c(0, 0),
                                breaks = c(if (up) h + br, if (down) -h - br), labels = c(if (up) br, if (down) br)) +
    theme_rodrigo(base_size) +
    ggplot2::theme(axis.title.x = ggplot2::element_blank(), axis.text.x = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(), panel.border = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_line(linewidth = 0.25),
                   legend.key.size = ggplot2::unit(3, "mm"))
}

#' Region of interest in detail (Taylor 2008 Fig. 4C)
#'
#' Summary score of every unified-breakpoint-profile region in a window with its analytical error, the score
#' thresholds for significance (`fdr`) and for peak detection (`peak_q`). The left axis gives q at the two
#' thresholds, halfway up the score range, and at the top bin of each peak. Underneath, each significant event is
#' a thin line (omitted when the event covers the window; its cytobands are then in the x title) and each peak a
#' line from boundary to boundary with a box on its top bin and one bin either side. Cancer genes in the window are
#' named along the top, each joined to the error-bar tip of the nearest region.
#'
#' @inheritParams rae_plot_genome
#' @param chrom one chromosome or arm, e.g. `13` or `"13q"`.
#' @param start,end window in bp; default the whole chromosome (or arm).
#' @param direction `"loss"` (D') or `"gain"` (A').
#' @return a ggplot.
#' @export
rae_plot_locus <- function(x, chrom, start = NULL, end = NULL, direction = c("loss", "gain"), genes = "oncokb",
                           focal_max_bp = 3e6, focal_max_genes = 10, germline_symbol = "Г", base_size = 6) {
  .need("ggplot2")
  direction <- match.arg(direction)
  L <- .layout(x, chrom[1])
  if (is.null(start)) start <- L$from
  if (is.null(end)) end <- L$to
  r <- x$regions
  key <- if (direction == "gain") "A" else "D"
  sc <- r[[paste0(key, "_score")]]; er <- r[[paste0(key, "_err")]]; q <- r[[paste0(key, "_q")]]
  w <- which(r$chrom == L$chrom & r$end >= start & r$start <= end)
  if (!length(w)) stop("no regions in the window")
  d <- data.frame(mid = (pmax(r$start[w], start) + pmin(r$end[w], end)) / 2 / 1e6, s = sc[w], e = er[w], q = q[w],
                  start = r$start[w], end = r$end[w])
  roi <- x$roi[x$roi$direction == direction & x$roi$chrom == L$chrom & x$roi$end >= start & x$roi$start <= end, , drop = FALSE]
  if (is.null(roi$germline_overlap)) roi$germline_overlap <- rep(FALSE, nrow(roi))
  pk <- roi[roi$type == "peak", , drop = FALSE]
  d$grp <- "region"
  for (i in seq_len(nrow(pk))) d$grp[d$start >= pk$start[i] & d$end <= pk$end[i]] <- if (pk$germline_overlap[i]) "germline" else direction
  qa <- .q_axis(sc, q, x$params$fdr)
  thr <- data.frame(y = c(qa$at[1], .q_axis(sc, q, x$params$peak_q)$at[1]),
                    lab = c(paste0("'FDR '*", x$params$fdr), paste0("'peak '*", x$params$peak_q)))
  thr <- thr[is.finite(thr$y), , drop = FALSE]
  lo <- min(d$s - d$e, thr$y, na.rm = TRUE); hi <- max(d$s + d$e, thr$y, na.rm = TRUE); rng <- hi - lo
  cb <- rae_cytobands(x$build)
  gt <- .genes_table(genes, x$build)
  d <- d[order(d$mid), , drop = FALSE]
  # events: a thin line each, unless the event covers the window (then it is named in the x title instead)
  ev <- roi[roi$type == "region", , drop = FALSE]
  covers <- ev$start <= start & ev$end >= end
  title_band <- if (any(covers)) .band_label(cb, L$chrom, ev$start[which(covers)[1]], ev$end[which(covers)[1]]) else .band_label(cb, L$chrom, start, end)
  ev <- ev[!covers, , drop = FALSE]
  ntr <- (nrow(ev) > 0) + (nrow(pk) > 0)
  tr <- function(k) lo - rng * (0.07 + 0.09 * (k - 1))   # track lines below the data
  tracks <- rbind(
    if (nrow(ev)) data.frame(type = "event", start = ev$start, end = ev$end, y = tr(1),
                             col = ifelse(ev$germline_overlap, "germline", "region"),
                             lab = ifelse(ev$germline_overlap, germline_symbol,
                                          vapply(seq_len(nrow(ev)), function(i) .band_label(cb, L$chrom, ev$start[i], ev$end[i]), character(1)))),
    if (nrow(pk)) data.frame(type = "peak", start = pk$start, end = pk$end, y = tr(ntr),
                             col = ifelse(pk$germline_overlap, "germline", direction),
                             lab = vapply(seq_len(nrow(pk)), function(i) {
                               if (pk$germline_overlap[i]) return(germline_symbol)
                               gsub("\n", ", ", .event_label(data.frame(chrom = L$chrom, start = pk$start[i], end = pk$end[i],
                                                                        size_bp = pk$end[i] - pk$start[i] + 1, germline_overlap = FALSE),
                                                             cb, gt, focal_max_bp, focal_max_genes, germline_symbol))
                             }, character(1))))
  pal <- c(rae_colours, region = rae_colours[["ink2"]])
  # peak: the most significant bin with one bin either side, as a box on the boundary-to-boundary line
  top <- NULL
  for (i in seq_len(nrow(pk))) {
    j <- which(d$start >= pk$start[i] & d$end <= pk$end[i])
    if (!length(j)) next
    k <- j[which.max(d$s[j])]
    top <- rbind(top, data.frame(xs = max(d$start[max(k - 1, 1)], start) / 1e6, xe = min(d$end[min(k + 1, nrow(d))], end) / 1e6,
                                 y = tr(ntr), s = d$s[k], q = d$q[k], col = tracks$col[tracks$type == "peak"][i]))
  }
  bottom <- if (ntr) tr(ntr) - rng * 0.08 else lo - rng * 0.03
  ytop <- hi + rng * (if (is.null(gt)) 0.03 else 0.22)
  p <- ggplot2::ggplot(d, ggplot2::aes(mid, s)) +
    ggplot2::geom_hline(data = thr, ggplot2::aes(yintercept = y), linetype = "22", linewidth = 0.25, colour = rae_colours[["ink2"]]) +
    ggplot2::geom_line(linewidth = 0.25, colour = rae_colours[["ink2"]]) +
    ggplot2::geom_errorbar(ggplot2::aes(ymin = s - e, ymax = s + e, colour = grp), width = 0, linewidth = 0.25) +
    ggplot2::geom_point(ggplot2::aes(colour = grp), size = 0.8) +
    ggplot2::scale_colour_manual(values = c(rae_colours, region = "#000000"), breaks = c(direction, "germline"),
                                 labels = c("peak", paste("overlaps germline variation", germline_symbol)), name = NULL)
  if (!is.null(tracks) && nrow(tracks)) {
    tracks$xs <- pmax(tracks$start, start) / 1e6; tracks$xe <- pmin(tracks$end, end) / 1e6
    p <- p + ggplot2::geom_segment(data = tracks, ggplot2::aes(x = xs, xend = xe, y = y, yend = y), colour = unname(pal[tracks$col]),
                                   linewidth = 0.4, lineend = "butt", inherit.aes = FALSE) +
      ggplot2::geom_text(data = tracks, ggplot2::aes(x = (xs + xe) / 2, y = y - rng * 0.025, label = paste(type, lab)), vjust = 1,
                         size = base_size / ggplot2::.pt, colour = rae_colours[["ink2"]], inherit.aes = FALSE)
  }
  if (!is.null(top))
    p <- p + ggplot2::geom_rect(data = top, ggplot2::aes(xmin = xs, xmax = xe, ymin = y - rng * 0.015, ymax = y + rng * 0.015),
                                fill = unname(pal[top$col]), colour = NA, inherit.aes = FALSE)
  # cancer genes: horizontal names along the top, each joined to the error-bar tip of the nearest region
  if (!is.null(gt)) {
    g <- gt[as.character(gt$chrom) == as.character(L$chrom) & gt$end >= start & gt$start <= end, , drop = FALSE]
    if (nrow(g) && nrow(g) <= 40) {
      g$x <- (pmax(g$start, start) + pmin(g$end, end)) / 2e6
      nn <- vapply(g$x, function(v) which.min(abs(d$mid - v)), integer(1))
      g$y <- d$s[nn] + d$e[nn]
      .need("ggrepel")
      p <- p + ggrepel::geom_text_repel(data = g, ggplot2::aes(x = x, y = y, label = gene), fontface = "italic",
                                        size = base_size / ggplot2::.pt, nudge_y = hi + rng * 0.14 - g$y, direction = "x",
                                        vjust = 0, segment.size = 0.15, segment.colour = rae_colours[["ns"]],
                                        min.segment.length = 0, box.padding = 0.1, max.overlaps = Inf, seed = 1,
                                        inherit.aes = FALSE)
    }
  }
  # left axis: the two thresholds, the q halfway up the visible score range, and the q at each peak's top bin
  q_at <- function(y) { k <- is.finite(sc) & is.finite(q) & sc >= y; if (any(k)) max(q[k]) else NA_real_ }
  fq <- function(v) { e <- floor(log10(v)); m <- signif(v / 10^e, 2)
    if (v >= 1e-3) paste0("'", signif(v, 2), "'") else paste0(m, "%*%10^", e) }
  ax <- thr
  mid <- (lo + hi) / 2
  if (is.finite(qm <- q_at(mid))) ax <- rbind(ax, data.frame(y = mid, lab = fq(qm)))
  if (!is.null(top)) ax <- rbind(ax, data.frame(y = top$s, lab = vapply(top$q, fq, character(1))))
  keep <- rep(TRUE, nrow(ax))
  for (i in seq_len(nrow(ax))[-1]) keep[i] <- all(abs(ax$y[i] - ax$y[seq_len(i - 1)][keep[seq_len(i - 1)]]) >= 0.05 * (ytop - bottom))
  ax <- ax[keep, , drop = FALSE]
  p + ggplot2::scale_y_continuous(name = paste0("summary score (", if (direction == "gain") "A'" else "D'", ")"),
                                  breaks = (function(b) b[b >= lo & b <= hi])(pretty(c(lo, hi))), position = "right",
                                  sec.axis = ggplot2::sec_axis(~ ., name = "q (FDR); dashed: thresholds",
                                                               breaks = ax$y, labels = parse(text = ax$lab))) +
    ggplot2::coord_cartesian(xlim = c(start, end) / 1e6, ylim = c(bottom, ytop), expand = FALSE) +
    ggplot2::labs(x = paste0("HSA ", L$chrom, " (Mb; ", title_band, ")")) +
    theme_rodrigo(base_size)
}

#' Copy-number heatmap of the cohort
#'
#' Wrapper around [ComplexHeatmap::Heatmap()] in the style of `HeatmapCNR`: genomic bins as rows in
#' chromosome order with the chromosome labelled on the left, samples as columns clustered by Bray-Curtis
#' dissimilarity (on 2^log2 ratio) with Ward's method (`ward.D2`). Values are the mode-centred log2 segment
#' means (`value = "log2"`) or the discrete calls (`value = "state"`, -2..2).
#'
#' @param x a `rae_result`.
#' @param chrom `NULL` (genome) or chromosomes/arms as in [rae_plot_genome()].
#' @param start,end optional window (bp) when a single chromosome is given (Taylor 2008 Fig. 4A).
#' @param value `"log2"` or `"state"`.
#' @param bin_size bin width in bp (default 1 Mb; 50 kb for windows).
#' @param col colour mapping (a `circlize::colorRamp2` function or named vector); default blue-white-red on
#'   log2 ratio -1..1 with white for |log2| < 0.16.
#' @param cluster_samples cluster columns (needs 'vegan').
#' @param base_size font size (points).
#' @param ... passed to [ComplexHeatmap::Heatmap()].
#' @return a `Heatmap`; draw it, or save with [rae_save()].
#' @export
rae_heatmap <- function(x, chrom = NULL, start = NULL, end = NULL, value = c("log2", "state"),
                        bin_size = if (is.null(start)) 1e6 else 5e4, col = NULL, cluster_samples = TRUE,
                        base_size = 6, ...) {
  .need(c("ComplexHeatmap", "circlize"))
  value <- match.arg(value)
  L <- .layout(x, chrom)
  if (!is.null(start) || !is.null(end)) {
    if (nrow(L) != 1) stop("start/end need a single chromosome")
    L$from <- max(L$from, if (is.null(start)) L$from else start); L$to <- min(L$to, if (is.null(end)) L$to else end)
  }
  bins <- do.call(rbind, lapply(seq_len(nrow(L)), function(k) {
    s <- seq(L$from[k], L$to[k], by = bin_size)
    data.frame(chrom = L$chrom[k], unit = L$label[k], start = s, end = pmin(s + bin_size - 1, L$to[k]))
  }))
  sg <- x$segments
  ids <- unique(x$samples$ID)
  v <- if (value == "log2") sg$seg.mean.centred else sg$state
  m <- matrix(NA_real_, nrow(bins), length(ids), dimnames = list(NULL, ids))
  mid <- (bins$start + bins$end) / 2
  for (j in seq_along(ids)) {
    s <- sg[sg$ID == ids[j], ]; vj <- v[sg$ID == ids[j]]
    for (c0 in unique(bins$chrom)) {
      i <- which(bins$chrom == c0); k <- which(s$chrom == c0)
      if (!length(k)) next
      o <- k[order(s$loc.start[k])]
      f <- findInterval(mid[i], s$loc.start[o])
      ok <- f > 0
      hit <- ok; hit[ok] <- mid[i][ok] <= s$loc.end[o][f[ok]]
      m[i[hit], j] <- vj[o][f[hit]]
    }
  }
  if (is.null(col)) col <- if (value == "log2")
    circlize::colorRamp2(c(-1, -0.16, 0.16, 1), c("#2166AC", "#FFFFFF", "#FFFFFF", "#D40000")) else
    c(`-2` = "#08306B", `-1` = "#4A90D9", `0` = "#FFFFFF", `1` = "#F4A3A3", `2` = "#D40000")
  if (value == "state") m[] <- as.character(m)
  cf <- factor(bins$unit, levels = L$label)
  mids <- vapply(split(seq_len(nrow(bins)), cf), function(i) i[ceiling(length(i) / 2)], numeric(1))
  chl <- stats::setNames(rep(c("#404040", "#BABABA"), length.out = nlevels(cf)), levels(cf))
  ra <- ComplexHeatmap::rowAnnotation(
    labs = ComplexHeatmap::anno_mark(at = mids, labels = names(mids), side = "left",
                                     labels_gp = grid::gpar(fontsize = base_size)),
    chr = cf, col = list(chr = chl), show_annotation_name = FALSE, show_legend = FALSE)
  cl <- FALSE
  if (cluster_samples && length(ids) > 2) {
    .need("vegan")
    mm <- if (value == "log2") 2^m else matrix(as.numeric(m) + 2, nrow(m))
    mm[is.na(mm)] <- if (value == "log2") 1 else 2
    cl <- stats::hclust(vegan::vegdist(t(mm), method = "bray"), method = "ward.D2")
  }
  ComplexHeatmap::Heatmap(m, name = if (value == "log2") "log2 ratio" else "call", col = col,
                          cluster_rows = FALSE, cluster_columns = cl, left_annotation = ra,
                          row_split = if (nrow(L) > 1) NULL else NULL, show_row_names = FALSE,
                          column_names_gp = grid::gpar(fontsize = base_size),
                          show_column_names = length(ids) <= 60, na_col = "#F0F0F0",
                          heatmap_legend_param = list(title_gp = grid::gpar(fontsize = base_size),
                                                      labels_gp = grid::gpar(fontsize = base_size)),
                          use_raster = TRUE, ...)
}

#' Noise estimate per tumour (derivative noise, DN)
#'
#' DN is the median absolute difference between neighbouring probes or bins (Taylor 2008 Eq. 1; the same
#' statistic as MAPD). It is measured when probe/bin values were supplied (`sample_stats`, or `bin_noise = TRUE`
#' for QDNAseq input) and otherwise estimated from the segmentation (gap G-01b, flag `dn_from_sd`).
#'
#' @param x a `rae_result`.
#' @param against `"rank"` (tumours ordered by DN) or `"E"` (the detector locations E of single-copy gain and
#'   loss against DN, showing how the per-tumour parameterisation follows the noise).
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_noise <- function(x, against = c("rank", "E"), base_size = 6) {
  .need("ggplot2")
  against <- match.arg(against)
  s <- x$samples
  s$source <- ifelse(grepl("dn_from_sd", s$flags), "estimated from segments", "measured from probes/bins")
  if (against == "rank") {
    s <- s[order(s$dn), ]; s$rank <- seq_len(nrow(s))
    return(ggplot2::ggplot(s, ggplot2::aes(rank, dn, shape = source)) +
      ggplot2::geom_hline(yintercept = stats::median(s$dn), linetype = "22", linewidth = 0.25, colour = rae_colours[["ink2"]]) +
      ggplot2::geom_point(size = 0.7) +
      ggplot2::scale_shape_manual(values = c(`measured from probes/bins` = 16, `estimated from segments` = 1), name = NULL) +
      ggplot2::scale_x_continuous(expand = c(0, 0)) +
      ggplot2::labs(x = "tumours, ordered by noise", y = "DN, median |difference| of neighbouring log2 ratios") +
      theme_rodrigo(base_size))
  }
  d <- rbind(data.frame(dn = s$dn, E = s$E_A0, k = "single-copy gain (E A0)", fb = grepl("A0_fallback", s$flags)),
             data.frame(dn = s$dn, E = s$E_D0, k = "hemizygous loss (E D0)", fb = grepl("D0_fallback", s$flags)))
  d$fb <- ifelse(d$fb, "fixed-slope fallback", "half-maximum")
  ggplot2::ggplot(d, ggplot2::aes(dn, E, colour = k, shape = fb)) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25, colour = rae_colours[["ns"]]) +
    ggplot2::geom_point(size = 0.7) +
    ggplot2::scale_colour_manual(values = c(`single-copy gain (E A0)` = rae_colours[["gain"]],
                                            `hemizygous loss (E D0)` = rae_colours[["loss"]]), name = NULL) +
    ggplot2::scale_shape_manual(values = c(`half-maximum` = 16, `fixed-slope fallback` = 4), name = NULL) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "DN (log2 ratio)", y = "detector location E (log2 ratio)") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2), shape = ggplot2::guide_legend(nrow = 2)) +
    theme_rodrigo(base_size)
}
