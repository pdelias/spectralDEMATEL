# The structure map: naming a system's type, and saying how firmly.
#
# A type is a verdict on a boundary, so it is worth exactly as much as the
# distance from that boundary. Everything here returns the margin alongside the
# name, because an interface cannot degrade its language without one.

#' The two cuts that partition the (coupling, hierarchy) plane.
#'
#' They are not the same kind of constant, and an interface must not present
#' them as though they were.
#'
#'   coupling  0.5 is the indirect-dominance threshold. At mu = 0.5 the
#'             dominant eigenvalue of the total-relation matrix is exactly 1, so
#'             indirect effects exactly equal direct ones. It follows from the
#'             algebra and is not a free parameter.
#'
#'   hierarchy 0.10 is a RECOMMENDATION. It sits near the corpus 75th
#'             percentile of 0.098 and reproduces the two-cluster solution on
#'             the two named types. A user may move it, and the interface should
#'             say that they may.
#'
#' @noRd
STRUCTURE_CUTS <- list(coupling = 0.5, hierarchy = 0.10)

#' Aggregate constants from the 117-system reference corpus.
#'
#' Aggregates only. The per-system rows are deliberately not shipped: they are
#' traceable to individual published studies and belong to work that is not yet
#' out. Everything the map and the trade-off reading need is here.
#' @noRd
CORPUS <- list(
  n = 117L,
  # Least-squares fit of hierarchy_sd on mu_max across the 117.
  tradeoff = list(intercept = 0.211766, slope = -0.204723,
                  resid_sd = 0.034985, r2 = 0.591, spearman = -0.776),
  # Counts by type, and the median total-to-direct multiplier within each.
  count = c(`diffuse-amplified` = 81L, `hierarchical-dampened` = 16L,
            `hierarchical-amplified` = 12L, `diffuse-dampened` = 8L),
  multiplier = c(`diffuse-amplified` = 5.14, `hierarchical-dampened` = 1.50,
                 `hierarchical-amplified` = 3.06, `diffuse-dampened` = 1.83),
  # 5th, 25th, 50th, 75th and 95th percentiles.
  quantiles = list(
    mu_max       = c(0.2432, 0.5824, 0.7541, 0.8555, 0.9187),
    hierarchy_sd = c(0.0140, 0.0279, 0.0549, 0.0985, 0.1686),
    dominance    = c(0.0121, 0.0253, 0.0457, 0.1032, 0.3389)),
  probs = c(0.05, 0.25, 0.50, 0.75, 0.95)
)

#' The intervention logic each type favours, in the source paper's own words.
#'
#' Quoted rather than paraphrased. Each is a hypothesis derived from structure,
#' not a validated result, and the wording was chosen to say exactly that much
#' and no more.
#' @noRd
INTERVENTION_LOGIC <- c(
  `diffuse-amplified` = paste(
    "A portfolio of coordinated changes, appraised on total effects.",
    "No factor offers materially more leverage than another."),
  `hierarchical-amplified` = paste(
    "Targeted action at the dominant factor, appraised on total",
    "effects, since most of what follows arrives indirectly."),
  `hierarchical-dampened` = paste(
    "Targeted action at the dominant factor, appraised on direct",
    "effects, which carry most of what follows."),
  `diffuse-dampened` = paste(
    "Effects stay near where they are applied. With neither",
    "multiplication nor a preferred point of entry, the system may",
    "divide into parts that can be treated separately.")
)

#' What the pairing of a type with an intervention logic is, and is not.
#' @noRd
INTERVENTION_CAVEAT <- paste(
  "Each pairing of a type with an intervention logic is a hypothesis about",
  "intervention design, derived from structure. Testing one takes outcome data",
  "of a kind the reference corpus does not hold, since it records the structure",
  "of expert causal models and not what followed from them. The hypothesis is",
  "specific -- it names the effects an option should be appraised on, and",
  "whether effort belongs on one factor or across many -- but it has not been",
  "validated against outcomes.")

#' Name a system's structural type, and say how firmly
#'
#' Places a diagnosed system on the (coupling, hierarchy) map and returns the
#' type, the margin to each cut in that axis's own units, and the intervention
#' logic the source paper pairs with the type.
#'
#' @param d A diagnostics list from [spectral_diagnostics()].
#' @param coupling_cut Where amplified meets dampened. The default 0.5 is the
#'   indirect-dominance threshold and follows from the algebra; changing it
#'   makes the type mean something other than what the corpus results mean.
#' @param hierarchy_cut Where hierarchical meets diffuse. The default 0.10 is a
#'   recommendation, not a fitted constant, and a user may move it.
#'
#' @return A list, or `NULL` when the diagnostics could not be computed:
#'   \describe{
#'     \item{`type`}{One of `"diffuse-amplified"`, `"hierarchical-amplified"`,
#'       `"hierarchical-dampened"`, `"diffuse-dampened"`.}
#'     \item{`coupling_margin`}{Signed distance from `coupling_cut`. Positive is
#'       amplified.}
#'     \item{`hierarchy_margin`}{Signed distance from `hierarchy_cut`. Positive
#'       is hierarchical.}
#'     \item{`nearest_cut`}{Which boundary the system sits closest to, once each
#'       margin is expressed in units of that axis's corpus spread.}
#'     \item{`nearest_margin_scaled`}{That distance, in corpus interquartile
#'       ranges. Small means the type is a close call.}
#'     \item{`corpus_share`}{Proportion of the 117 reference systems of this
#'       type, so a user in an uncommon corner is told so.}
#'     \item{`corpus_multiplier`}{Median total-to-direct multiplier within the
#'       type across the corpus.}
#'     \item{`intervention_logic`}{The source paper's wording for this type.}
#'     \item{`caveat`}{What that pairing is and is not. Show it wherever the
#'       logic is shown.}
#'     \item{`cuts`}{The two cuts actually used, and whether each was the
#'       default.}
#'     \item{`hierarchy_reading`}{Always `"hierarchy_sd"`. See below.}
#'   }
#'
#' @section Only the standard-deviation reading classifies:
#' The 0.10 cut is calibrated against `hierarchy_sd`, and every corpus-level
#' result is expressed in it. `hierarchy_gini` and `hierarchy_pr` measure the
#' same idea on different scales -- and the participation ratio runs the
#' opposite way -- so a cut of 0.10 means nothing on either. Show all three
#' readings by all means; classify on one.
#'
#' @section The margin is the point:
#' A system at coupling 0.51 and one at 0.95 are both amplified, and the
#' confidence in the advice is not remotely the same. This function returns the
#' distance and leaves the wording to the caller, because how much to hedge is a
#' presentation decision.
#'
#' @section The advice is a hypothesis:
#' `intervention_logic` is quoted from the source paper, where each pairing is
#' explicitly a hypothesis derived from structure and not a validated result.
#' An interface that renders it as an instruction will be read as stronger than
#' the evidence. `caveat` exists to be shown, not stored.
#'
#' @examples
#' d <- spectral_diagnostics(worked_matrices$fefo_stock_control)
#' s <- structural_type(d)
#' s$type
#' s$coupling_margin
#' s$intervention_logic
#'
#' # An uncommon corner is worth knowing about:
#' structural_type(spectral_diagnostics(worked_matrices$resilience_capabilities))$corpus_share
#'
#' @seealso [type_stability()] for whether the type survives moving the cut,
#'   [tradeoff_residual()] for how the system sits against the corpus trade-off.
#' @export
structural_type <- function(d,
                            coupling_cut  = STRUCTURE_CUTS$coupling,
                            hierarchy_cut = STRUCTURE_CUTS$hierarchy) {
  if (is.null(d) || !is.finite(d$mu_max) || !is.finite(d$hierarchy_sd)) return(NULL)

  type <- paste0(if (d$hierarchy_sd > hierarchy_cut) "hierarchical" else "diffuse",
                 "-",
                 if (d$mu_max > coupling_cut) "amplified" else "dampened")

  coupling_margin  <- d$mu_max - coupling_cut
  hierarchy_margin <- d$hierarchy_sd - hierarchy_cut

  # The two axes are not on comparable scales -- coupling spans (0,1) and
  # hierarchy_sd rarely exceeds 0.2 -- so "which cut is this system closest to"
  # is meaningless until both distances are expressed in the same currency. The
  # corpus interquartile range of each axis is that currency.
  iqr <- function(q) q[4] - q[2]
  scaled <- c(coupling  = abs(coupling_margin)  / iqr(CORPUS$quantiles$mu_max),
              hierarchy = abs(hierarchy_margin) / iqr(CORPUS$quantiles$hierarchy_sd))
  nearest <- names(scaled)[which.min(scaled)]

  list(
    type                  = type,
    coupling_margin       = coupling_margin,
    hierarchy_margin      = hierarchy_margin,
    nearest_cut           = nearest,
    nearest_margin_scaled = unname(scaled[nearest]),
    corpus_share          = unname(CORPUS$count[type]) / CORPUS$n,
    corpus_multiplier     = unname(CORPUS$multiplier[type]),
    intervention_logic    = unname(INTERVENTION_LOGIC[type]),
    caveat                = INTERVENTION_CAVEAT,
    cuts = list(coupling  = coupling_cut,
                hierarchy = hierarchy_cut,
                coupling_is_default  = identical(coupling_cut,  STRUCTURE_CUTS$coupling),
                hierarchy_is_default = identical(hierarchy_cut, STRUCTURE_CUTS$hierarchy)),
    hierarchy_reading = "hierarchy_sd"
  )
}

#' Does the type survive moving the hierarchy cut?
#'
#' The hierarchy cut is a recommendation rather than a fitted constant, so the
#' honest question is not "which side is this system on" but "does the answer
#' change across the range anyone might reasonably choose". For one matrix that
#' is a handful of classifications and a sentence.
#'
#' @param d A diagnostics list from [spectral_diagnostics()].
#' @param hierarchy_cuts Cuts to try. The default spans the range the source
#'   paper's supplement re-runs its whole corpus assignment across.
#' @param coupling_cut Held fixed: it is not a free parameter.
#'
#' @return A list, or `NULL` when the diagnostics could not be computed:
#'   \describe{
#'     \item{`by_cut`}{Data frame of `cut` and the resulting `type`.}
#'     \item{`stable`}{`TRUE` when every cut gives the same type.}
#'     \item{`flips_between`}{The two adjacent cuts the type changes between,
#'       or `NA` when it does not change.}
#'   }
#'
#' @examples
#' d <- spectral_diagnostics(worked_matrices$fefo_stock_control)
#' s <- type_stability(d)
#' s$stable
#' s$by_cut
#'
#' @export
type_stability <- function(d,
                           hierarchy_cuts = c(0.075, 0.09, 0.10, 0.11),
                           coupling_cut = STRUCTURE_CUTS$coupling) {
  if (is.null(d) || !is.finite(d$hierarchy_sd)) return(NULL)

  hierarchy_cuts <- sort(unique(hierarchy_cuts))
  types <- vapply(hierarchy_cuts, function(h) {
    structural_type(d, coupling_cut = coupling_cut, hierarchy_cut = h)$type
  }, character(1))

  changed <- which(types[-1] != types[-length(types)])
  list(
    by_cut = data.frame(cut = hierarchy_cuts, type = types,
                        stringsAsFactors = FALSE),
    stable = length(changed) == 0L,
    flips_between = if (length(changed) == 0L) NA_real_ else
      c(hierarchy_cuts[changed[1]], hierarchy_cuts[changed[1] + 1L])
  )
}

#' Where the system sits against the corpus coupling-hierarchy trade-off
#'
#' Coupling and hierarchy oppose each other across the reference corpus,
#' coupling accounting for about 59% of the variation in hierarchy. A system's
#' residual from that relation says whether it has unusually concentrated
#' leverage for how strongly it is coupled.
#'
#' @param d A diagnostics list from [spectral_diagnostics()].
#'
#' @return A list, or `NULL` when the diagnostics could not be computed:
#'   `expected` (hierarchy predicted from this system's coupling), `residual`,
#'   `residual_sd` (that residual in units of the corpus residual standard
#'   deviation), `direction` (`"above"` or `"below"` the line), and `caveat`.
#'
#' @section A supplement, never a substitute:
#' Used on its own to assign types, the trade-off is wrong for 19 of the 117
#' reference systems, and wrong where it counts: 10 systems with a dominant
#' point of entry would be told to spread their effort, and all 8 without one
#' would be sent to look for one. Read this beside the computed hierarchy, never
#' instead of it.
#'
#' @examples
#' d <- spectral_diagnostics(worked_matrices$resilience_capabilities)
#' r <- tradeoff_residual(d)
#' round(r$residual_sd, 2)
#' r$direction
#'
#' @export
tradeoff_residual <- function(d) {
  if (is.null(d) || !is.finite(d$mu_max) || !is.finite(d$hierarchy_sd)) return(NULL)

  tr <- CORPUS$tradeoff
  expected <- tr$intercept + tr$slope * d$mu_max
  residual <- d$hierarchy_sd - expected

  list(
    expected    = expected,
    residual    = residual,
    residual_sd = residual / tr$resid_sd,
    direction   = if (residual >= 0) "above" else "below",
    caveat = paste(
      "The trade-off states what to expect of a system with this much",
      "coupling; the computed hierarchy states what this system actually has.",
      "Used alone to assign types it is wrong for 19 of the 117 reference",
      "systems, and wrong where it counts. Read it beside the computed",
      "hierarchy, never instead of it.")
  )
}
