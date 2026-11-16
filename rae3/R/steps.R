#' Individual steps of RAE, exposed for validation
#'
#' Thin wrappers around single functions of the C++ core, so each step of the method can be checked
#' on its own rather than only end to end.
#'
#' * `rae_sigmoid()`: the soft discriminator, Eq. 2 of the paper. For `component = "A1"` the
#'   truncated, rescaled form is used.
#' * `rae_dsigmoid()`: its first derivative as used in the analytical error (Methods S1-C).
#' * `rae_density()`: Gaussian kernel density on `n` points over `[from, to]`, reproducing
#'   `stats::density()`; `rae_bw_nrd0()` reproduces `stats::bw.nrd0()`.
#' * `rae_bh()`: Benjamini-Hochberg q-values, identical to `p.adjust(method = "BH")`.
#' * `rae_empirical_p()`: `(b + 1) / (N + 1)` with `b` the number of null scores at or above the
#'   observed score.
#' * `rae_rng_unif()`: uniform draws from the core's own generator.
#' * `rae_builds()`: names of the bundled genome builds.
#'
#' @param a numeric vector of log2 ratios.
#' @param E,beta location and slope of the detector.
#' @param component one of `"A0"`, `"A1"`, `"D0"`, `"D1"`.
#' @param x numeric vector.
#' @param bw bandwidth (standard deviation of the kernel).
#' @param from,to,n range and number of evaluation points.
#' @param p numeric vector of p-values.
#' @param observed observed scores.
#' @param null_scores pooled null scores.
#' @param seed,stream seed and sub-stream of the generator.
#' @return Numeric vectors, except `rae_density()` (list with `x`, `y`, `bw`) and `rae_builds()`
#'   (character).
#' @name rae_steps
#' @export rae_sigmoid rae_dsigmoid rae_density rae_bw_nrd0 rae_bh rae_empirical_p rae_rng_unif rae_builds
#' @aliases rae_sigmoid rae_dsigmoid rae_density rae_bw_nrd0 rae_bh rae_empirical_p rae_rng_unif rae_builds
#' @usage
#' rae_sigmoid(a, E, beta, component = "A0")
#' rae_dsigmoid(a, E, beta, component = "A0")
#' rae_density(x, bw, from = -1, to = 1, n = 2048L)
#' rae_bw_nrd0(x)
#' rae_bh(p)
#' rae_empirical_p(observed, null_scores)
#' rae_rng_unif(n, seed, stream = 0)
#' rae_builds()
NULL
