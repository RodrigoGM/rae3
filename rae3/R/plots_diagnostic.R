# Diagnostic figures of Taylor et al. 2008 (PLoS ONE e3179), recomputed from a rae_result.
#   Fig. 2C, S3, S5 (right)  rae_plot_density()        per-tumour density and the four detectors
#   Fig. S5 (left)           rae_plot_sample()         per-tumour genome profile
#   Fig. S10                 rae_plot_E()              detector locations across the cohort
#   Fig. S9                  rae_plot_diploid_offset() detected diploid offset vs median of the segmentation
#   Fig. S4                  rae_plot_entropy()        segment count vs entropy change (needs probe/bin data)
#   Fig. 3A, 4B              rae_plot_segments()       per-tumour segments, UBP boundaries, hotspots in a window
#   Fig. 3B                  rae_plot_hotspots()       original vs randomised hotspot density on one arm
#   Fig. 3C                  rae_plot_region_sizes()   region sizes, UBP vs hotspot-cleaved
#   (not in the paper)       rae_plot_arms()           arm-level test (G-14e): score vs arm length, per-arm bars
# Not reproduced: Fig. 1 and 2A-B (schematics), S2 (HapMap genotypes), S6-S8 (data-specific examples, LOH).

.sample_index <- function(x, sample) {
  i <- if (is.numeric(sample)) sample else match(sample, x$samples$ID)
  if (is.na(i) || i < 1 || i > nrow(x$samples)) stop("sample not found: ", sample)
  i
}

#' Per-tumour density and detectors (Taylor 2008 Fig. 2C, S3, S5)
#'
#' The marker-weighted density of the tumour's segment means (scaled to a maximum of 1), with the four soft
#' detectors: single-copy gain A0, amplification A1, hemizygous loss D0 and homozygous deletion D1, and their
#' fitted location E and slope beta. The density is recomputed with the run's parameters and random stream, so
#' it is the one the parameters were fitted on.
#'
#' @param x a `rae_result`.
#' @param sample sample ID or row number in `x$samples`.
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_density <- function(x, sample, base_size = 6) {
  .need("ggplot2")
  i <- .sample_index(x, sample)
  s <- x$samples[i, ]
  sg <- x$segments[x$segments$ID == s$ID, ]
  fit <- rae_parameterize(sg$seg.mean.centred + s$shift, sg$seg.sd, sg$num.mark, dn = s$dn, params = x$params, stream = i - 1)
  dx <- fit$density_x - fit$shift
  d <- data.frame(x = dx, y = fit$density_y / max(fit$density_y))
  grid <- seq(-1, 1, length.out = 801)
  comp <- c("A0", "A1", "D0", "D1")
  lab <- c(A0 = "single-copy gain (A0)", A1 = "amplification (A1)", D0 = "hemizygous loss (D0)", D1 = "homozygous deletion (D1)")
  sig <- do.call(rbind, lapply(comp, function(k)
    data.frame(x = grid, y = rae_sigmoid(grid, s[[paste0("E_", k)]], s[[paste0("beta_", k)]], k), k = lab[[k]])))
  sig$k <- factor(sig$k, levels = lab)
  txt <- paste(sprintf("%s  E = %.3f  β = %.1f", comp, unlist(s[paste0("E_", comp)]), unlist(s[paste0("beta_", comp)])), collapse = "\n")
  ggplot2::ggplot() +
    ggplot2::geom_vline(xintercept = 0, linewidth = 0.25, colour = rae_colours[["ns"]]) +
    ggplot2::geom_line(data = d, ggplot2::aes(x, y), linewidth = 0.35, colour = rae_colours[["ink"]]) +
    ggplot2::geom_line(data = sig, ggplot2::aes(x, y, colour = k, linetype = k), linewidth = 0.35) +
    ggplot2::scale_colour_manual(values = stats::setNames(rae_colours[c("gain", "gain", "loss", "loss")], lab), name = NULL) +
    ggplot2::scale_linetype_manual(values = stats::setNames(c("solid", "42", "solid", "12"), lab), name = NULL) +
    ggplot2::labs(x = "log2 ratio (mode-centred)", y = "density (scaled to 1) / detector value", caption = txt,
                  title = paste0(s$ID, "   DN = ", signif(s$dn, 3), if (nzchar(s$flags) && s$flags != "-") paste0("   flags: ", s$flags))) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(-1, 1), ylim = c(0, 1.02)) +
    theme_rodrigo(base_size) + ggplot2::theme(plot.title = ggplot2::element_text(size = base_size),
                                            plot.caption = ggplot2::element_text(hjust = 0, size = base_size, lineheight = 0.9))
}

.bin_values <- function(bins) {
  if (is.null(bins)) return(NULL)
  if (inherits(bins, "QDNAseqCopyNumbers")) {
    fd <- Biobase::fData(bins); cn <- Biobase::assayDataElement(bins, "copynumber")
    ch <- sub("^chr", "", as.character(fd$chromosome)); keep <- ch %in% as.character(1:22) & (if (is.null(fd$use)) TRUE else fd$use)
    return(data.frame(chrom = as.integer(ch[keep]), start = fd$start[keep], end = fd$end[keep],
                      log2(pmax(cn[keep, , drop = FALSE], 0) + .Machine$double.xmin), check.names = FALSE))
  }
  if (!is.data.frame(bins) || ncol(bins) < 4) stop("bins must be a QDNAseqCopyNumbers object or a data frame chrom, start, end, <log2 ratio per sample>")
  names(bins)[1:3] <- c("chrom", "start", "end"); bins$chrom <- as.integer(sub("^chr", "", bins$chrom))
  bins[bins$chrom %in% 1:22, ]
}

#' Per-tumour genome profile (Taylor 2008 Fig. S5, left)
#'
#' Mode-centred log2 segment means across the genome, coloured by the discrete call (gain, loss, neutral);
#' with `bins`, the probe/bin log2 ratios are drawn underneath (shifted by the same centring).
#'
#' @inheritParams rae_plot_density
#' @param chrom `NULL` or chromosomes/arms, as in [rae_plot_genome()].
#' @param bins optional probe/bin values: a `QDNAseqCopyNumbers` object or a data frame chrom, start, end and
#'   one log2-ratio column per sample.
#' @return a ggplot.
#' @export
rae_plot_sample <- function(x, sample, chrom = NULL, bins = NULL, base_size = 6) {
  .need("ggplot2")
  i <- .sample_index(x, sample); s <- x$samples[i, ]
  L <- .layout(x, chrom)
  sg <- x$segments[x$segments$ID == s$ID, ]
  iv <- .ginterval(L, sg$chrom, sg$loc.start, sg$loc.end)
  seg <- data.frame(xmin = iv$xmin, xmax = iv$xmax, y = sg$seg.mean.centred[iv$row],
                    call = ifelse(sg$state[iv$row] > 0, "gain", ifelse(sg$state[iv$row] < 0, "loss", "neutral")))
  p <- ggplot2::ggplot() +
    ggplot2::geom_rect(data = data.frame(xmin = L$offset, xmax = L$offset + L$to - L$from)[c(FALSE, TRUE), , drop = FALSE],
                       ggplot2::aes(xmin = xmin, xmax = xmax, ymin = -Inf, ymax = Inf), fill = "#F2F2F2")
  bv <- .bin_values(bins)
  if (!is.null(bv)) {
    if (!s$ID %in% names(bv)) stop("sample ", s$ID, " not among the bin columns")
    b <- data.frame(x = .gpos(L, bv$chrom, (bv$start + bv$end) / 2), y = bv[[s$ID]] - s$shift)
    b <- b[is.finite(b$x) & is.finite(b$y), ]
    p <- p + ggplot2::geom_point(data = b, ggplot2::aes(x, y), size = 0.05, colour = rae_colours[["ns"]], stroke = 0)
  }
  p + ggplot2::geom_hline(yintercept = 0, linewidth = 0.2, colour = rae_colours[["ink2"]]) +
    ggplot2::geom_segment(data = seg, ggplot2::aes(x = xmin, xend = xmax, y = y, yend = y, colour = call), linewidth = 0.6) +
    ggplot2::scale_colour_manual(values = c(gain = rae_colours[["gain"]], loss = rae_colours[["loss"]], neutral = "#000000"), name = NULL) +
    ggplot2::scale_x_continuous(breaks = L$offset + (L$to - L$from) / 2, labels = L$label, expand = c(0, 0)) +
    ggplot2::coord_cartesian(ylim = c(-1.5, 1.5)) +
    ggplot2::labs(x = NULL, y = "log2 ratio (mode-centred)", title = s$ID) +
    theme_rodrigo(base_size) + ggplot2::theme(plot.title = ggplot2::element_text(size = base_size))
}

#' Detector locations across the cohort (Taylor 2008 Fig. S10)
#'
#' E of single-copy gain (A0, bars up) and hemizygous loss (D0, bars down) for every tumour, each sorted;
#' tumours whose detector fell back to the fixed slope are marked.
#'
#' @param x a `rae_result`.
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_E <- function(x, base_size = 6) {
  .need("ggplot2")
  s <- x$samples; n <- nrow(s)
  a <- s[order(s$E_A0), ]; d <- s[order(-s$E_D0), ]
  df <- rbind(data.frame(i = seq_len(n), E = a$E_A0, k = "gain", fb = grepl("A0_fallback", a$flags)),
              data.frame(i = n + seq_len(n), E = d$E_D0, k = "loss", fb = grepl("D0_fallback", d$flags)))
  ggplot2::ggplot(df, ggplot2::aes(i, E, fill = k)) +
    ggplot2::geom_col(width = 1) +
    ggplot2::geom_point(data = df[df$fb, ], ggplot2::aes(i, E * 1.04), shape = 4, size = 0.5, stroke = 0.3, show.legend = FALSE) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = c(gain = rae_colours[["gain"]], loss = rae_colours[["loss"]]),
                               labels = c(gain = "single-copy gain, E A0", loss = "hemizygous loss, E D0"), name = NULL) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = paste0("tumours (n = ", n, "), each component sorted;  × fixed-slope fallback"),
                  y = "detector location E (log2 ratio)") +
    theme_rodrigo(base_size) + ggplot2::theme(axis.text.x = ggplot2::element_blank(), axis.ticks.x = ggplot2::element_blank())
}

#' Detected diploid offset (Taylor 2008 Fig. S9)
#'
#' The centring shift found from each tumour's density (the diploid peak) against the marker-weighted median
#' of its original segment means. Points off the identity line are tumours whose diploid peak is not the
#' median, e.g. after whole-genome doubling or with a large altered fraction.
#'
#' @param x a `rae_result`.
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_diploid_offset <- function(x, base_size = 6) {
  .need("ggplot2")
  s <- x$samples; sg <- x$segments
  med <- vapply(seq_len(nrow(s)), function(i) {
    g <- sg[sg$ID == s$ID[i], ]; v <- g$seg.mean.centred + s$shift[i]
    o <- order(v); w <- cumsum(g$num.mark[o]) / sum(g$num.mark)
    v[o][which(w >= 0.5)[1]]
  }, numeric(1))
  ggplot2::ggplot(data.frame(m = med, sh = s$shift), ggplot2::aes(m, sh)) +
    ggplot2::geom_abline(slope = 1, intercept = 0, linetype = "22", linewidth = 0.25, colour = rae_colours[["gain"]]) +
    ggplot2::geom_point(size = 0.7) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "median log2 ratio of the segmentation (marker-weighted)", y = "detected diploid peak offset (log2 ratio)") +
    theme_rodrigo(base_size)
}

#' Segment count against entropy change (Taylor 2008 Fig. S4)
#'
#' Entropy change Delta S = S(r) - S(m) from the probe/bin log2 ratios r to the segmented values m (each probe
#' carrying its segment's mean), with S the Shannon entropy of a 2048-point density on -1..1 (as stated in the
#' paper). The paper's bandwidth is "fixed" but not given; `bw` is therefore a choice (gap register G-02c).
#' Hyper-segmented tumours have many segments and a small Delta S.
#'
#' @param x a `rae_result`.
#' @param bins probe/bin values (required), as in [rae_plot_sample()].
#' @param bw density bandwidth (log2 units). Not stated in the paper.
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_entropy <- function(x, bins, bw = 0.02, base_size = 6) {
  .need("ggplot2")
  bv <- .bin_values(bins)
  H <- function(v) { d <- stats::density(v[is.finite(v)], bw = bw, from = -1, to = 1, n = 2048)$y; p <- d / sum(d); -sum(p[p > 0] * log(p[p > 0])) }
  s <- x$samples
  res <- do.call(rbind, lapply(seq_len(nrow(s)), function(i) {
    id <- s$ID[i]; if (!id %in% names(bv)) return(NULL)
    sg <- x$segments[x$segments$ID == id, ]
    r <- bv[[id]]; m <- rep(NA_real_, length(r))
    for (k in seq_len(nrow(sg))) {
      j <- bv$chrom == sg$chrom[k] & bv$start >= sg$loc.start[k] & bv$end <= sg$loc.end[k]
      m[j] <- sg$seg.mean.centred[k] + s$shift[i]
    }
    data.frame(ID = id, dS = H(r) - H(m), n = nrow(sg))
  }))
  if (is.null(res)) stop("no sample IDs match the bin columns")
  ggplot2::ggplot(res, ggplot2::aes(dS, n)) +
    ggplot2::geom_smooth(method = "lm", formula = y ~ x, linewidth = 0.3, colour = rae_colours[["gain"]], fill = "#F4C7C7") +
    ggplot2::geom_point(size = 0.7, shape = 1) +
    ggplot2::scale_y_log10() +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = "Δ entropy, probes/bins to segments", y = "number of segments") +
    theme_rodrigo(base_size)
}

#' Segments of each tumour in a window (Taylor 2008 Fig. 3A and 4B)
#'
#' One line per tumour carrying a call in the window: thin for single-copy changes (gain or hemizygous loss),
#' thick for high-level ones (amplification or homozygous deletion). Underneath, the boundaries of the
#' unified breakpoint profile; on top, recombination-hotspot centres used by the null model.
#'
#' @inheritParams rae_plot_locus
#' @return a ggplot.
#' @export
rae_plot_segments <- function(x, chrom, start = NULL, end = NULL, direction = c("loss", "gain"), base_size = 6) {
  .need("ggplot2")
  direction <- match.arg(direction)
  L <- .layout(x, chrom[1])
  if (is.null(start)) start <- L$from
  if (is.null(end)) end <- L$to
  sg <- x$segments
  sel <- sg$chrom == L$chrom & sg$loc.end >= start & sg$loc.start <= end & (if (direction == "loss") sg$state < 0 else sg$state > 0)
  g <- sg[sel, ]
  ids <- unique(g$ID)
  g$row <- match(g$ID, ids)
  g$level <- ifelse(abs(g$state) >= 2, if (direction == "loss") "homozygous deletion" else "amplification",
                    if (direction == "loss") "hemizygous loss" else "single-copy gain")
  r <- x$regions[x$regions$chrom == L$chrom & x$regions$start >= start & x$regions$start <= end, ]
  hs <- x$diagnostics$hotspots
  hs <- if (is.null(hs)) NULL else hs[hs$chrom == L$chrom & hs$pos >= start & hs$pos <= end, ]
  n <- max(1, length(ids))
  p <- ggplot2::ggplot() +
    ggplot2::geom_segment(data = g, ggplot2::aes(x = pmax(loc.start, start) / 1e6, xend = pmin(loc.end, end) / 1e6, y = row, yend = row,
                                                 linewidth = level), colour = .dir_col(direction), lineend = "butt") +
    ggplot2::scale_linewidth_manual(values = stats::setNames(c(0.3, 1.0), if (direction == "loss") c("hemizygous loss", "homozygous deletion") else c("single-copy gain", "amplification")), name = NULL) +
    ggplot2::geom_segment(data = r, ggplot2::aes(x = start / 1e6, xend = start / 1e6, y = -0.5, yend = -0.5 - 0.04 * n), linewidth = 0.15)
  if (!is.null(hs) && nrow(hs)) p <- p + ggplot2::geom_segment(data = hs, ggplot2::aes(x = pos / 1e6, xend = pos / 1e6, y = n + 1, yend = n + 1 + 0.04 * n),
                                                               linewidth = 0.15, colour = rae_colours[["ink2"]])
  p + ggplot2::annotate("text", x = start / 1e6, y = c(-0.5 - 0.06 * n, n + 1 + 0.06 * n),
                        label = c("unified breakpoint profile (region starts)", if (!is.null(hs) && nrow(hs)) "recombination hotspots" else ""),
                        hjust = 0, vjust = c(1, 0), size = base_size / ggplot2::.pt, colour = rae_colours[["ink2"]]) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = c(start, end) / 1e6, ylim = c(-0.5 - 0.12 * n, n + 1 + 0.12 * n)) +
    ggplot2::labs(x = paste0("chromosome ", L$label, " (Mb)"), y = paste0("tumours with a call (n = ", length(ids), ")")) +
    theme_rodrigo(base_size) + ggplot2::theme(axis.text.y = ggplot2::element_blank(), axis.ticks.y = ggplot2::element_blank())
}

#' Original and randomised recombination hotspots (Taylor 2008 Fig. 3B)
#'
#' Density of the hotspot centres on one chromosome arm and of the cut positions drawn from them in the first
#' randomisation of the null model (rejection sampling, spec step 8). The two should have the same shape.
#'
#' @param x a `rae_result`.
#' @param chrom one chromosome or arm, e.g. `"1p"`.
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_hotspots <- function(x, chrom = "1p", base_size = 6) {
  .need("ggplot2")
  d <- x$diagnostics
  if (is.null(d) || !isTRUE(d$has_hotspots)) stop("this result has no hotspot map (custom build or older result)")
  L <- .layout(x, chrom[1])
  pick <- function(t) t$pos[t$chrom == L$chrom & t$pos >= L$from & t$pos <= L$to]
  h <- pick(d$hotspots); cuts <- pick(d$sampled_cuts)
  bw <- stats::median(diff(sort(h)))
  dens <- function(v, k) { e <- stats::density(v, bw = bw, from = L$from, to = L$to, n = 1024); data.frame(x = e$x / 1e6, y = e$y / max(e$y), k = k) }
  df <- rbind(dens(h, paste0("hotspot centres (n = ", length(h), ")")), dens(cuts, paste0("randomised cuts, first round (n = ", length(cuts), ")")))
  rug <- rbind(data.frame(x = h / 1e6, y = -0.04, k = df$k[1]), data.frame(x = cuts / 1e6, y = -0.10, k = utils::tail(df$k, 1)))
  ggplot2::ggplot(df, ggplot2::aes(x, y, colour = k)) +
    ggplot2::geom_line(linewidth = 0.3) +
    ggplot2::geom_point(data = rug, shape = "|", size = 0.8) +
    ggplot2::scale_colour_manual(values = stats::setNames(c(rae_colours[["ink"]], rae_colours[["loss"]]), unique(df$k)), name = NULL) +
    ggplot2::scale_x_continuous(expand = c(0, 0)) +
    ggplot2::labs(x = paste0("chromosome ", L$label, " (Mb)"), y = paste0("density (scaled to 1); bandwidth ", round(bw / 1e3), " kb")) +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    theme_rodrigo(base_size)
}

#' Region sizes before and after hotspot cleavage (Taylor 2008 Fig. 3C)
#'
#' @param x a `rae_result`.
#' @param base_size font size (points).
#' @return a ggplot.
#' @export
rae_plot_region_sizes <- function(x, base_size = 6) {
  .need("ggplot2")
  d <- x$diagnostics
  if (is.null(d$ubp_bp)) stop("this result has no null-model diagnostics")
  df <- rbind(data.frame(bp = d$ubp_bp, k = "unified breakpoint profile"),
              if (length(d$cleaved_bp)) data.frame(bp = d$cleaved_bp, k = "hotspot-cleaved (first null round)"))
  df <- df[df$bp > 0, , drop = FALSE]
  # 50 log-spaced bins over the data; both axes end exactly at the data (largest region, tallest bin)
  rx <- range(log10(df$bp)); ed <- seq(rx[1], rx[2], length.out = 51)
  h <- do.call(rbind, lapply(split(df, df$k), function(z) {
    n <- tabulate(findInterval(log10(z$bp), ed, rightmost.closed = TRUE, all.inside = TRUE), 50)
    data.frame(xmin = 10^ed[-51], xmax = 10^ed[-1], n = n, k = z$k[1])
  }))
  ymax <- max(h$n)
  yb <- pretty(c(0, ymax)); yb <- c(yb[yb <= ymax * 0.92], ymax)
  xmax <- 10^rx[2]
  xb <- 10^(ceiling(rx[1]):floor(rx[2])); xb <- c(xb[log10(xb) <= rx[2] - 0.3], xmax)
  size_lab <- function(v) ifelse(v >= 1e6, paste(signif(v / 1e6, 2), "Mb"), ifelse(v >= 1e3, paste(signif(v / 1e3, 2), "kb"), paste(signif(v, 2), "bp")))
  ggplot2::ggplot(h, ggplot2::aes(xmin = xmin, xmax = xmax, ymin = 0, ymax = n, colour = k, fill = k)) +
    ggplot2::geom_rect(linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = c(`unified breakpoint profile` = "#8F8F8F", `hotspot-cleaved (first null round)` = NA), name = NULL) +
    ggplot2::scale_colour_manual(values = c(`unified breakpoint profile` = "#404040", `hotspot-cleaved (first null round)` = rae_colours[["loss"]]), name = NULL) +
    ggplot2::scale_x_log10(breaks = xb, labels = size_lab(xb), expand = c(0, 0)) +
    ggplot2::scale_y_continuous(breaks = yb, expand = c(0, 0)) +
    ggplot2::coord_cartesian(xlim = 10^rx, ylim = c(0, ymax)) +
    ggplot2::labs(x = "size of region", y = "number of regions") +
    ggplot2::guides(colour = ggplot2::guide_legend(nrow = 2)) +
    theme_rodrigo(base_size) +
    ggplot2::theme(plot.margin = ggplot2::margin(3, 6, 1, 1, unit = "mm"))   # room for the end labels
}

#' Arm-level cohort test (gap register G-14e)
#'
#' The table of the arm-level test (`arm_test = TRUE` in [rae_params()]). `type = "length"` plots each arm's
#' cohort gain score A′ = sqrt(A0² + A1²) and loss score D′ = sqrt(D0² + D1²) against the arm's length in Mb:
#' the test treats arms as exchangeable, so a trend with length warns that longer arms have a different
#' background rate. `type = "bar"` draws one bar per arm in genome order, gain up and loss down; arms with q at
#' or below `fdr` are coloured, the rest grey.
#'
#' @param x a `rae_result` from a run with `arm_test = TRUE`.
#' @param type `"length"` or `"bar"`.
#' @param base_size font size (points).
#' @param fdr q-value cut-off for the highlight; by default the run's `fdr`.
#' @param ... unused.
#' @return a ggplot.
#' @export
rae_plot_arms <- function(x, type = c("length", "bar"), base_size = 6, fdr = x$params$fdr, ...) {
  .need("ggplot2")
  type <- match.arg(type)
  a <- x$arms
  if (is.null(a)) stop("this result has no arm-level table: run with arm_test = TRUE in rae_params()", call. = FALSE)
  lab <- c(gain = "gain (A′)", loss = "loss (D′)")
  long <- function(k) data.frame(arm = a$arm_label, i = seq_len(nrow(a)), len = a$length_bp / 1e6,
                                 score = a[[paste0("score_", k)]], q = a[[paste0("q_", k)]], k = k)
  df <- rbind(long("gain"), long("loss"))
  df$sig <- is.finite(df$q) & df$q <= fdr
  cols <- c(gain = rae_colours[["gain"]], loss = rae_colours[["loss"]])
  if (type == "length") {
    return(ggplot2::ggplot(df, ggplot2::aes(len, score, colour = k)) +
      ggplot2::geom_point(ggplot2::aes(shape = sig), size = 0.9, stroke = 0.3) +
      ggplot2::geom_text(data = df[df$sig, , drop = FALSE], ggplot2::aes(label = arm), size = base_size / ggplot2::.pt,
                         vjust = -0.7, show.legend = FALSE) +
      ggplot2::scale_colour_manual(values = cols, labels = lab, name = NULL) +
      ggplot2::scale_shape_manual(values = c(`FALSE` = 1, `TRUE` = 16),
                                  labels = c(`FALSE` = paste("q >", fdr), `TRUE` = paste("q ≤", fdr)), name = NULL, drop = FALSE) +
      ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.03, 0.08))) +   # room for the labels
      ggplot2::labs(x = "arm length (Mb)", y = "cohort arm score") +
      theme_rodrigo(base_size))
  }
  df$y <- ifelse(df$k == "gain", df$score, -df$score)
  df$fill <- ifelse(df$sig, df$k, "ns")
  ggplot2::ggplot(df, ggplot2::aes(i, y, fill = fill)) +
    ggplot2::geom_col(width = 0.8) +
    ggplot2::geom_hline(yintercept = 0, linewidth = 0.25) +
    ggplot2::scale_fill_manual(values = c(cols, ns = rae_colours[["ns"]]),
                               labels = c(gain = paste("gain, q ≤", fdr), loss = paste("loss, q ≤", fdr), ns = "not significant"),
                               name = NULL) +
    ggplot2::scale_x_continuous(breaks = seq_len(nrow(a)), labels = a$arm_label, expand = c(0.01, 0)) +
    ggplot2::scale_y_continuous(labels = abs) +
    ggplot2::labs(x = NULL, y = "cohort arm score (gain up, loss down)") +
    theme_rodrigo(base_size) +
    ggplot2::theme(axis.text.x = ggplot2::element_text(angle = 90, vjust = 0.5, hjust = 1))
}
