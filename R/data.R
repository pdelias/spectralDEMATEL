#' Two transcribed DEMATEL influence matrices
#'
#' The two worked systems from the source paper's Section 5, transcribed from
#' their published tables and verified against those papers' own printed
#' total-relation matrices. They are the package's regression fixtures and the
#' worked examples in the vignette.
#'
#' Between them they straddle every threshold the diagnostics use: coupling
#' above and below 0.5, hierarchy above and below 0.10, indirect dominance both
#' ways, and strong connectivity both ways.
#'
#' @format A named list of two numeric matrices. Each carries a `doi` attribute
#'   and a `label` attribute.
#'
#'   \describe{
#'     \item{`fefo_stock_control`}{13 factors, DOI `10.1016/j.jafr.2025.101848`.
#'       Strongly coupled (`mu_max` 0.828), diffuse, single-mode. Has a
#'       **non-zero diagonal**, which makes it a useful case for input-form
#'       detection.}
#'     \item{`resilience_capabilities`}{15 factors, DOI
#'       `10.1007/s12063-024-00470-8`. Weakly coupled (`mu_max` 0.351),
#'       hierarchical. **Fails strong connectivity**: row 14 is all zeros, so
#'       factor 14 dispatches nothing and assumption A2 does not hold. The
#'       source paper reports its diagnostics regardless, which is why an
#'       assumption check has to be a returned verdict rather than a refusal.}
#'   }
#'
#' @source Transcribed in the JDS paper repository, `data/raw/matrices/`.
#'
#' @examples
#' names(worked_matrices)
#' attr(worked_matrices$fefo_stock_control, "doi")
#'
#' d <- spectral_diagnostics(worked_matrices$fefo_stock_control)
#' d$mu_max
#'
#' # The other one is not strongly connected:
#' is_irreducible(worked_matrices$resilience_capabilities)
"worked_matrices"
