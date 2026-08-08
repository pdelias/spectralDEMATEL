# How much a structural type is worth.
#
# A type is a verdict on a boundary. structural_type() says how far the system
# sits from that boundary; these two functions ask the other questions. Would
# the same reading arise from the rating distribution alone? And does it survive
# the noise expert ratings are known to carry?

#' Where the observed diagnostics sit in their own surrogate ensemble
#'
#' Shuffles the off-diagonal entries, holding the number of factors, the density
#' and the exact multiset of values fixed, and asks where the real system falls
#' in the resulting null. A diagnosis is worth less if the same numbers would
#' arise from the rating distribution on its own.
#'
#' This is the comparison the source paper could not run across its corpus,
#' because the precomputed spectral file holds no raw matrices. A user of this
#' package has their own matrix, so they can.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#' @param B Ensemble size. 200 is the paper's figure; lower it when a person is
#'   waiting.
#' @param seed Optional integer, so a figure in a paper can be reproduced.
#' @param coupling_cut,hierarchy_cut Passed to [structural_type()] when
#'   classifying each draw.
#'
#' @return A list, or `NULL` when the matrix cannot be diagnosed:
#'   \describe{
#'     \item{`metrics`}{Data frame, one row per diagnostic, with the observed
#'       value, the ensemble median, minimum and maximum, `share_ge` (the
#'       proportion of draws at or above the observed value) and `outside`
#'       (`TRUE` when the observation falls beyond the whole ensemble).}
#'     \item{`type`}{The observed structural type.}
#'     \item{`type_share`}{Proportion of surrogate draws that produce the same
#'       type. High means the type follows from the rating distribution rather
#'       than from the structure the analyst built.}
#'     \item{`type_table`}{Counts of every type across the ensemble.}
#'     \item{`B`}{Ensemble size actually used.}
#'     \item{`requested`}{Ensemble size asked for.}
#'     \item{`complete`}{`TRUE` when `B` reached `requested`.}
#'     \item{`attempts`}{Shuffles drawn, accepted and rejected together.}
#'   }
#'
#'   Returns `NULL` when the matrix is inadmissible, when it is not irreducible,
#'   or when no shuffle of it was admissible.
#'
#' @section Read `complete` before quoting a share:
#' A sparse matrix rejects most shuffles, so the ensemble can come back shorter
#' than requested; see [surrogate_ensemble()]. A `share_ge` computed over seven
#' draws and one computed over two hundred are not the same claim, and `B` alone
#' does not distinguish them from what the caller asked for.
#'
#' @section Reading `share_ge`:
#' It is a position, not a p-value. `share_ge = 0` means the observation exceeds
#' every draw; `share_ge = 1` means every draw exceeds it. Values near 0.5 mean
#' the shuffle reproduces the observation, so that diagnostic says more about
#' how the ratings were distributed than about how the factors were connected.
#'
#' @section `type_share` is weak evidence for a common type:
#' Nearly 70% of the reference corpus is diffuse-amplified, so a shuffled matrix
#' usually lands there too. A `type_share` near 1 for a common type therefore
#' means very little -- it says the null also produces that type, not that the
#' observation is unremarkable. The per-metric positions are the informative
#' part. `type_share` earns its keep when the observed type is an uncommon one,
#' where a high value is a genuine warning that the rating distribution alone
#' would have put the system in that corner.
#'
#' @examples
#' sp <- surrogate_position(worked_matrices$fefo_stock_control, B = 40, seed = 1)
#' sp$metrics[, c("metric", "observed", "share_ge", "outside")]
#' sp$type_share
#'
#' @param max_attempts Maximum number of shuffles drawn before giving up,
#'   default `50 * B`, passed to [surrogate_ensemble()]. Reaching it yields a
#'   short ensemble and `complete = FALSE` rather than an error.
#'
#' @seealso [surrogate_ensemble()] for the raw draws.
#' @export
surrogate_position <- function(A, B = 200, seed = NULL,
                               coupling_cut  = STRUCTURE_CUTS$coupling,
                               hierarchy_cut = STRUCTURE_CUTS$hierarchy,
                               max_attempts  = 50L * B) {
  obs <- spectral_diagnostics(A, checks = FALSE)
  if (is.null(obs)) return(NULL)
  if (!is_irreducible(A)) return(NULL)

  ens <- surrogate_ensemble(A, B = B, seed = seed, max_attempts = max_attempts)
  # No admissible shuffle at all. There is no null distribution to sit in, so
  # there is nothing to report -- same answer shape as an inadmissible matrix.
  if (is.null(ens)) return(NULL)

  quantities <- c("mu_max", "hierarchy_sd", "dominance")
  metrics <- do.call(rbind, lapply(quantities, function(q) {
    e <- ens[[q]]
    data.frame(
      metric   = q,
      observed = obs[[q]],
      median   = stats::median(e),
      min      = min(e),
      max      = max(e),
      share_ge = mean(e >= obs[[q]]),
      outside  = obs[[q]] < min(e) || obs[[q]] > max(e),
      stringsAsFactors = FALSE)
  }))
  rownames(metrics) <- NULL

  obs_type <- structural_type(obs, coupling_cut, hierarchy_cut)$type
  draw_types <- vapply(seq_len(nrow(ens)), function(i) {
    structural_type(list(mu_max = ens$mu_max[i], hierarchy_sd = ens$hierarchy_sd[i]),
                    coupling_cut, hierarchy_cut)$type
  }, character(1))

  list(
    metrics    = metrics,
    type       = obs_type,
    type_share = mean(draw_types == obs_type),
    type_table = table(draw_types),
    B          = nrow(ens),
    # What was asked for, versus what the rejection rate allowed. A share
    # computed over seven draws and one computed over two hundred are not the
    # same claim, and nothing downstream can tell them apart from `B` alone.
    requested  = B,
    complete   = isTRUE(attr(ens, "complete")),
    attempts   = attr(ens, "attempts")
  )
}

#' Does the type survive the noise expert ratings carry?
#'
#' Perturbs the influence ratings by a stated tolerance and recomputes the
#' structural type. Expert judgements on a four- or five-point scale carry at
#' least half a point of noise, and a classification that does not survive that
#' is not a classification.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#' @param tolerance Half-width of the perturbation, in rating-scale points.
#'   The default 0.5 is the noise a four- or five-point expert scale is
#'   generally taken to carry.
#' @param B Number of perturbed matrices to draw.
#' @param seed Optional integer, so a result can be reproduced.
#' @param perturb Which entries to disturb. `"nonzero"`, the default, leaves
#'   recorded zeros alone; `"all"` disturbs every off-diagonal entry.
#' @param coupling_cut,hierarchy_cut Passed to [structural_type()].
#'
#' @return A list, or `NULL` when the matrix cannot be diagnosed:
#'   `observed_type`, `share_same` (proportion of perturbed matrices keeping
#'   that type), `type_table`, `B_admissible` (draws that could be diagnosed at
#'   all), `tolerance` and `perturb`.
#'
#' @section Why recorded zeros are left alone by default:
#' A rating of zero is a judgement that there is no influence between two
#' factors, which is a different kind of statement from a rating of one carrying
#' measurement error. Disturbing zeros upward also fills in the matrix, which
#' destroys the sparsity the connectivity checks are about. `perturb = "all"` is
#' available for anyone who disagrees, and the choice is reported back.
#'
#' Perturbed values are clamped to `[0, max(A)]`: a rating cannot go negative,
#' and nothing should exceed the top of the scale the analyst used.
#'
#' @examples
#' ms <- measurement_stability(worked_matrices$fefo_stock_control,
#'                             B = 40, seed = 1)
#' ms$observed_type
#' ms$share_same
#'
#' @export
measurement_stability <- function(A, tolerance = 0.5, B = 200, seed = NULL,
                                  perturb = c("nonzero", "all"),
                                  coupling_cut  = STRUCTURE_CUTS$coupling,
                                  hierarchy_cut = STRUCTURE_CUTS$hierarchy) {
  perturb <- match.arg(perturb)
  obs <- spectral_diagnostics(A, checks = FALSE)
  if (is.null(obs)) return(NULL)
  observed_type <- structural_type(obs, coupling_cut, hierarchy_cut)$type

  if (!is.null(seed)) set.seed(seed)

  n <- nrow(A)
  off <- which(row(A) != col(A))
  targets <- if (perturb == "all") off else off[A[off] > 0]
  ceiling_value <- max(A)

  types <- character(0)
  for (b in seq_len(B)) {
    P <- A
    noise <- stats::runif(length(targets), -tolerance, tolerance)
    P[targets] <- pmin(pmax(P[targets] + noise, 0), ceiling_value)
    d <- spectral_diagnostics(P, checks = FALSE)
    if (is.null(d)) next
    st <- structural_type(d, coupling_cut, hierarchy_cut)
    if (is.null(st)) next
    types <- c(types, st$type)
  }

  list(
    observed_type = observed_type,
    share_same    = if (length(types)) mean(types == observed_type) else NA_real_,
    type_table    = table(types),
    B_admissible  = length(types),
    tolerance     = tolerance,
    perturb       = perturb
  )
}
