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
#'   hierarchy 0.10 is a RECOMMENDATION under both readings, for different
#'             reasons. On H it is a tolerance: under the decision model of the
#'             source paper, uniform effort forgoes the share H of the
#'             attainable value, so 0.10 is the loss accepted before a model
#'             counts as hierarchical. On sd(u) it sat near the upper quartile
#'             of the 117-system corpus. A user may move it, and the interface
#'             should say that they may.
#'
#' The two readings share the number 0.10 and not its meaning. For a unit-length
#' entry profile H = (n - 1) sd^2, so sd = 0.10 is H = 0.10 only at n = 11.
#'
#' @noRd
STRUCTURE_CUTS <- list(coupling = 0.5, hierarchy = 0.10)

#' The band around the H cut that the source paper reports its sensitivity on
#' (htype_bands in the paper's R/pipeline/config.R).
#' @noRd
HIERARCHY_BANDS <- c(0.05, 0.20)

#' The cuts that type_stability() tries by default under each reading.
#'
#' Under H: the ends of the band, the cut, and a quarter either side of it, as
#' hier_axis() in the paper's R/pipeline/config.R computes them. Under sd: the
#' grid of 0.6.0, which is the one the first version of the paper used.
#' @noRd
hierarchy_grid <- function(reading) {
  if (reading == "sd") return(c(0.075, 0.09, 0.10, 0.11))
  b <- sort(HIERARCHY_BANDS)
  sort(unique(round(c(b, STRUCTURE_CUTS$hierarchy * c(0.75, 1, 1.25)), 4)))
}

#' The diagnostics field each reading classifies on.
#' @noRd
reading_field <- function(reading) {
  switch(reading, H = "hierarchy_h", sd = "hierarchy_sd")
}

#' Both classifying readings run the same way. Returned with every type so
#' that the direction never has to be remembered.
#' @noRd
HIERARCHY_DIRECTION <- "high means concentrated"

#' One finite number.
#' @noRd
is_number <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x)
}

#' Aggregate constants of the two reference corpora, one per reading.
#'
#' Aggregates only. The per-system rows are deliberately not shipped: they are
#' traceable to individual published studies. Everything the map and the
#' trade-off reading need is here.
#'
#'   sd  the 117 systems of the first version of the source paper, classified
#'       on sd(u). Frozen at the values of 0.6.0, so that an earlier diagnosis
#'       can be reproduced.
#'   H   the 81 models of the revised paper, classified on H. Built by
#'       data-raw/corpus_constants.R from the paper's pipeline and held in
#'       R/corpus-h.R.
#'
#' The two corpora differ in more than the reading. The revised paper screened
#' the studies again, zeroes the diagonal of every matrix before normalising,
#' and fits its trade-off line by median regression.
#' @noRd
CORPUS <- list(
  sd = list(
    n = 117L,
    # Least-squares fit of hierarchy_sd on mu_max across the 117.
    tradeoff = list(intercept = 0.211766, slope = -0.204723,
                    resid_sd = 0.034985, r2 = 0.591, spearman = -0.776,
                    method = "least squares"),
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
  ),
  H = CORPUS_H
)

#' The intervention logic each type is paired with, in the source paper's words.
#'
#' Quoted, one set per reading. `sd` is the wording of the first version of the
#' paper, as 0.6.0 shipped it. `H` is the wording of quadrant_profile() in the
#' revised paper's R/results_access.R, with the symbol H* written out as "the
#' tolerance H*". Each is a hypothesis derived from structure and not a
#' validated result, and the wording was chosen to say that much.
#' @noRd
INTERVENTION_LOGIC <- list(
  sd = c(
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
  ),
  H = c(
    `diffuse-amplified` = paste(
      "Spread effort, appraised on total effects. Spreading forgoes at",
      "most the tolerance H* of the attainable value, direct appraisal",
      "more than a quarter."),
    `hierarchical-amplified` = paste(
      "Effort targeted by entry weight, appraised on total effects.",
      "Spreading forgoes more than the tolerance H*, direct appraisal",
      "more than a quarter."),
    `hierarchical-dampened` = paste(
      "Effort targeted by entry weight, appraised on direct effects.",
      "Spreading forgoes more than the tolerance H*, direct appraisal at",
      "most a quarter."),
    `diffuse-dampened` = paste(
      "Spread effort, appraised on direct effects. Spreading forgoes at",
      "most the tolerance H*, direct appraisal at most a quarter.")
  )
)

#' What the pairing of a type with an intervention logic is, and is not.
#'
#' `sd` is the caveat of 0.6.0. `H` follows the closing paragraph of Section 5.1
#' of the revised paper.
#' @noRd
INTERVENTION_CAVEAT <- list(
  sd = paste(
    "Each pairing of a type with an intervention logic is a hypothesis about",
    "intervention design, derived from structure. Testing one takes outcome data",
    "of a kind the reference corpus does not hold, since it records the structure",
    "of expert causal models and not what followed from them. The hypothesis is",
    "specific -- it names the effects an option should be appraised on, and",
    "whether effort belongs on one factor or across many -- but it has not been",
    "validated against outcomes."),
  H = paste(
    "Each pairing of a type with an intervention hypothesis follows from the",
    "two cuts under the assumptions of a decision model, and holds by",
    "construction inside that model. Those assumptions are restrictive. Unequal",
    "costs, a fixed budget or set-up costs change the share that uniform effort",
    "captures. Complementarities between interventions fall outside the",
    "additive model. Whether a pairing holds for the system a model describes",
    "is an empirical question. A test takes outcome data of a kind the",
    "reference corpus does not hold, since the corpus records the structure of",
    "expert causal models and not what followed from them. The hypothesis is",
    "specific. A type names the effects an appraisal would count and whether",
    "effort would follow the entry weights or spread across the factors. It",
    "has not been validated against outcomes.")
)

#' The assumptions behind the two effort shares.
#'
#' From Appendix B of the revised source paper (the model and Proposition B.1,
#' items iii and iv).
#' @noRd
EFFORT_ASSUMPTIONS <- paste(
  "Both shares hold inside one decision model. Effort on a factor has a",
  "quadratic cost. The value of the resulting change is linear, there is",
  "neither a budget nor a set-up cost, and the objective is the dominant mode",
  "of propagation. With equal costs across the factors, uniform effort forgoes",
  "the share H of the value attainable along the dominant mode. Appraisal on",
  "direct effects alone forgoes the share mu_max^2 whatever the costs, one",
  "quarter at a coupling of one half. Unequal costs change the first share. The",
  "model does not cover complementarities between interventions.")

#' Name a system's structural type, and say how firmly
#'
#' Places a diagnosed system on the (coupling, hierarchy) map and returns the
#' type, the margin to each cut in that axis's own units, and the intervention
#' logic the source paper pairs with the type.
#'
#' @param d A diagnostics list from [spectral_diagnostics()].
#' @param coupling_cut Where amplified meets dampened. The default 0.5 is the
#'   indirect-dominance threshold and follows from the algebra. Changing it
#'   makes the type mean something other than what the corpus results mean.
#' @param hierarchy_cut Where hierarchical meets diffuse, on the scale of
#'   `reading`. The default 0.10 is a recommendation under both readings and a
#'   user may move it.
#' @param reading Which hierarchy reading classifies. `"H"`, the default, uses
#'   `hierarchy_h`. `"sd"` uses `hierarchy_sd` and reproduces the types of
#'   version 0.6.0 and earlier. **Both run the same way: high means
#'   concentrated.**
#'
#' @return A list, or `NULL` when the diagnostics could not be computed, when
#'   they do not hold the field the reading needs, or when a cut is not a
#'   finite number:
#'   \describe{
#'     \item{`type`}{One of `"diffuse-amplified"`, `"hierarchical-amplified"`,
#'       `"hierarchical-dampened"`, `"diffuse-dampened"`.}
#'     \item{`coupling_margin`}{Signed distance from `coupling_cut`. Positive is
#'       amplified.}
#'     \item{`hierarchy_margin`}{Signed distance from `hierarchy_cut`, in the
#'       units of the reading. Positive is hierarchical.}
#'     \item{`nearest_cut`}{Which boundary the system sits closest to, once each
#'       margin is expressed in units of that axis's corpus spread.}
#'     \item{`nearest_margin_scaled`}{That distance, in corpus interquartile
#'       ranges. Small means the type is a close call.}
#'     \item{`corpus_share`}{Proportion of the reference corpus of this type,
#'       so a user in an uncommon corner is told so.}
#'     \item{`corpus_multiplier`}{Median total-to-direct multiplier within the
#'       type across the corpus.}
#'     \item{`intervention_logic`}{The source paper's wording for this type.}
#'     \item{`caveat`}{What that pairing is and is not. Show it wherever the
#'       logic is shown.}
#'     \item{`cuts`}{The two cuts actually used, and whether each was the
#'       default.}
#'     \item{`hierarchy_reading`}{`"hierarchy_h"` or `"hierarchy_sd"`: the
#'       field that was classified on.}
#'     \item{`corpus_n`}{Size of the reference corpus behind `corpus_share`:
#'       81 models under `"H"`, 117 systems under `"sd"`.}
#'     \item{`hierarchy_direction`}{`"high means concentrated"`, which holds
#'       for both readings.}
#'     \item{`effort`}{A list of `uniform_forgoes`, `direct_appraisal_forgoes`
#'       and `assumptions`. See below. Both shares are `NA` under
#'       `reading = "sd"`.}
#'   }
#'
#' @section Two readings classify, and they are not interchangeable:
#' `hierarchy_h` and `hierarchy_sd` are tied by \eqn{H = (n - 1)\,
#' \mathrm{sd}(u)^2}, so a cut of 0.10 on one is a cut of 0.10 on the other only
#' for eleven factors. A larger model can be hierarchical on H and diffuse on
#' sd, and a smaller one the reverse. The type returned by 0.7.0 can therefore
#' differ from the type 0.6.0 returned for the same matrix. To reproduce a
#' 0.6.0 type, call `structural_type(d, reading = "sd")`.
#'
#' Each reading is read against its own reference corpus: 81 models for H and
#' 117 systems for sd. `hierarchy_gini` and `hierarchy_pr` do not classify.
#' They measure the same idea on other scales, and the participation ratio runs
#' the opposite way.
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
#' The wording is written for the default cuts. Under `reading = "H"` it speaks
#' of "the tolerance H*", which is `hierarchy_cut`, and of "a quarter", which is
#' the share forgone at a coupling of one half. With another `coupling_cut` the
#' quarter no longer marks the boundary of the type, and
#' `cuts$coupling_is_default` is `FALSE` to say so.
#'
#' @section The effort shares:
#' Under `reading = "H"` the two cuts have a reading in units of value, inside
#' the decision model of the source paper. `uniform_forgoes` is the share of the
#' value attainable along the dominant mode that effort spread evenly across
#' the factors gives up, and it equals `hierarchy_h`. `direct_appraisal_forgoes`
#' is the share given up by appraising an option on its direct effects alone,
#' and it equals `mu_max^2`, one quarter at the coupling cut. Neither is a
#' forecast. `assumptions` states the model they hold in, and belongs wherever
#' either share is shown.
#'
#' @examples
#' d <- spectral_diagnostics(worked_matrices$fefo_stock_control)
#' s <- structural_type(d)
#' s$type
#' s$hierarchy_reading
#' s$coupling_margin
#' s$intervention_logic
#' s$effort$uniform_forgoes
#'
#' # The type that version 0.6.0 returned:
#' structural_type(d, reading = "sd")$type
#'
#' # An uncommon corner is worth knowing about:
#' structural_type(spectral_diagnostics(worked_matrices$resilience_capabilities))$corpus_share
#'
#' @seealso [type_stability()] for whether the type survives moving the cut,
#'   [tradeoff_residual()] for how the system sits against the corpus trade-off.
#' @export
structural_type <- function(d,
                            coupling_cut  = STRUCTURE_CUTS$coupling,
                            hierarchy_cut = STRUCTURE_CUTS$hierarchy,
                            reading = c("H", "sd")) {
  reading <- match.arg(reading)
  field <- reading_field(reading)
  if (!is.list(d) || !is_number(d$mu_max) || !is_number(d[[field]])) return(NULL)
  if (!is_number(coupling_cut) || !is_number(hierarchy_cut)) return(NULL)

  hierarchy <- d[[field]]
  cp <- CORPUS[[reading]]

  type <- paste0(if (hierarchy > hierarchy_cut) "hierarchical" else "diffuse",
                 "-",
                 if (d$mu_max > coupling_cut) "amplified" else "dampened")

  coupling_margin  <- d$mu_max - coupling_cut
  hierarchy_margin <- hierarchy - hierarchy_cut

  # The two axes are not on comparable scales, so "which cut is this system
  # closest to" is meaningless until both distances are expressed in the same
  # currency. The corpus interquartile range of each axis is that currency.
  iqr <- function(q) q[4] - q[2]
  scaled <- c(coupling  = abs(coupling_margin)  / iqr(cp$quantiles$mu_max),
              hierarchy = abs(hierarchy_margin) / iqr(cp$quantiles[[field]]))
  nearest <- names(scaled)[which.min(scaled)]

  list(
    type                  = type,
    coupling_margin       = coupling_margin,
    hierarchy_margin      = hierarchy_margin,
    nearest_cut           = nearest,
    nearest_margin_scaled = unname(scaled[nearest]),
    corpus_share          = unname(cp$count[type]) / cp$n,
    corpus_multiplier     = unname(cp$multiplier[type]),
    intervention_logic    = unname(INTERVENTION_LOGIC[[reading]][type]),
    caveat                = INTERVENTION_CAVEAT[[reading]],
    cuts = list(coupling  = coupling_cut,
                hierarchy = hierarchy_cut,
                coupling_is_default  = identical(coupling_cut,  STRUCTURE_CUTS$coupling),
                hierarchy_is_default = identical(hierarchy_cut, STRUCTURE_CUTS$hierarchy)),
    hierarchy_reading   = field,
    corpus_n            = cp$n,
    hierarchy_direction = HIERARCHY_DIRECTION,
    # Proposition B.1 of the source paper, items iii and iv. The shares are
    # stated on H, so they are not reported when the type was read on sd.
    effort = list(
      uniform_forgoes          = if (reading == "H") hierarchy else NA_real_,
      direct_appraisal_forgoes = if (reading == "H") d$mu_max^2 else NA_real_,
      assumptions              = EFFORT_ASSUMPTIONS)
  )
}

#' Does the type survive moving the hierarchy cut?
#'
#' The hierarchy cut is a recommendation, so the question that matters is
#' whether the type changes across the range a user might reasonably choose.
#' For one matrix that is a handful of classifications and a sentence.
#'
#' @param d A diagnostics list from [spectral_diagnostics()].
#' @param hierarchy_cuts Cuts to try, on the scale of `reading`. `NULL`, the
#'   default, uses the grid of the reading: 0.05, 0.075, 0.10, 0.125 and 0.20
#'   under `"H"` (the ends of the band the source paper reports, the cut, and a
#'   quarter either side of it), and 0.075, 0.09, 0.10 and 0.11 under `"sd"`
#'   (the grid of version 0.6.0).
#' @param coupling_cut Held fixed: it is not a free parameter.
#' @param reading Which hierarchy reading classifies, as in [structural_type()].
#'   **Both run the same way: high means concentrated.**
#'
#' @return A list, or `NULL` when the diagnostics could not be computed or the
#'   cuts are not finite numbers:
#'   \describe{
#'     \item{`by_cut`}{Data frame of `cut` and the resulting `type`.}
#'     \item{`stable`}{`TRUE` when every cut gives the same type.}
#'     \item{`flips_between`}{The two adjacent cuts the type changes between,
#'       or `NA` when it does not change.}
#'     \item{`hierarchy_reading`}{`"hierarchy_h"` or `"hierarchy_sd"`.}
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
                           hierarchy_cuts = NULL,
                           coupling_cut = STRUCTURE_CUTS$coupling,
                           reading = c("H", "sd")) {
  reading <- match.arg(reading)
  field <- reading_field(reading)
  if (!is.list(d) || !is_number(d$mu_max) || !is_number(d[[field]])) return(NULL)
  if (!is_number(coupling_cut)) return(NULL)
  if (is.null(hierarchy_cuts)) hierarchy_cuts <- hierarchy_grid(reading)
  if (!is.numeric(hierarchy_cuts) || length(hierarchy_cuts) == 0L ||
      !all(is.finite(hierarchy_cuts))) return(NULL)

  hierarchy_cuts <- sort(unique(hierarchy_cuts))
  types <- vapply(hierarchy_cuts, function(h) {
    structural_type(d, coupling_cut = coupling_cut, hierarchy_cut = h,
                    reading = reading)$type
  }, character(1))

  changed <- which(types[-1] != types[-length(types)])
  list(
    by_cut = data.frame(cut = hierarchy_cuts, type = types,
                        stringsAsFactors = FALSE),
    stable = length(changed) == 0L,
    flips_between = if (length(changed) == 0L) NA_real_ else
      c(hierarchy_cuts[changed[1]], hierarchy_cuts[changed[1] + 1L]),
    hierarchy_reading = field
  )
}

#' Where the system sits against the corpus coupling-hierarchy trade-off
#'
#' Coupling and hierarchy oppose each other across the reference corpus. A
#' system's residual from that relation says whether it has unusually
#' concentrated influence for how strongly it is coupled.
#'
#' @param d A diagnostics list from [spectral_diagnostics()].
#' @param reading Which hierarchy reading, and with it which corpus line: the
#'   line of the 81-model corpus on `hierarchy_h` under `"H"`, the default, or
#'   the line of the 117-system corpus on `hierarchy_sd` under `"sd"`. **Both
#'   readings run the same way: high means concentrated.**
#'
#' @return A list, or `NULL` when the diagnostics could not be computed:
#'   `expected` (hierarchy predicted from this system's coupling), `residual`,
#'   `residual_sd` (that residual in units of the corpus residual standard
#'   deviation), `direction` (`"above"` or `"below"` the line), `caveat`,
#'   `hierarchy_reading` and `fit` (how the corpus line was fitted).
#'
#' @section The two lines:
#' Under `"H"` the line is the median regression of the revised source paper,
#' \eqn{H = 0.374 - 0.436\,\mu_{max}}, and coupling accounts for about 75% of
#' the variation in H. The line is straight and H is not negative, so above a
#' coupling of about 0.86 `expected` falls below zero. Read the residual there
#' as a distance from the line and not from a typical model. Under `"sd"` the
#' line is the least-squares fit of version 0.6.0, with coupling accounting for
#' about 59% of the variation.
#'
#' @section A supplement, never a substitute:
#' Used on its own to assign types, the trade-off is wrong for 9 of the 81
#' reference models on H: 3 hierarchical-amplified models would be told to
#' spread their effort, and all 3 diffuse-dampened models would be classed as
#' hierarchical. On sd it is wrong for 19 of the 117 reference systems. Read
#' this beside the computed hierarchy, never instead of it.
#'
#' @examples
#' d <- spectral_diagnostics(worked_matrices$resilience_capabilities)
#' r <- tradeoff_residual(d)
#' round(r$residual_sd, 2)
#' r$direction
#'
#' @export
tradeoff_residual <- function(d, reading = c("H", "sd")) {
  reading <- match.arg(reading)
  field <- reading_field(reading)
  if (!is.list(d) || !is_number(d$mu_max) || !is_number(d[[field]])) return(NULL)

  cp <- CORPUS[[reading]]
  tr <- cp$tradeoff
  expected <- tr$intercept + tr$slope * d$mu_max
  residual <- d[[field]] - expected

  list(
    expected    = expected,
    residual    = residual,
    residual_sd = residual / tr$resid_sd,
    direction   = if (residual >= 0) "above" else "below",
    caveat = if (reading == "sd") paste(
      "The trade-off states what to expect of a system with this much",
      "coupling; the computed hierarchy states what this system actually has.",
      "Used alone to assign types it is wrong for 19 of the 117 reference",
      "systems, and wrong where it counts. Read it beside the computed",
      "hierarchy, never instead of it.") else paste(
      "The trade-off states what to expect of a model with this much coupling.",
      "The computed hierarchy states what this model has. Used alone to",
      "assign types, the trade-off is wrong for", cp$alone$wrong, "of the",
      cp$n, "reference models, and wrong where it counts. Read it beside the",
      "computed hierarchy, never instead of it."),
    hierarchy_reading = field,
    fit = tr$method
  )
}
