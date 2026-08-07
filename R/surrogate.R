# The permutation surrogate baseline.
#
# Ported verbatim from R/spectral_dematel.R in the JDS paper repository.

#' Permutation surrogate baseline
#'
#' Shuffles the off-diagonal entries of `A` while holding the number of factors,
#' the density and the exact multiset of entry values fixed, then recomputes the
#' diagnostics. Repeating this gives a null distribution against which an
#' observed diagnostic can be read.
#'
#' The question it answers is the one a single number cannot: how much of the
#' structure is something the analyst built, and how much follows from the
#' entry distribution on its own. A system whose coupling sits in the middle of
#' its own surrogate ensemble has coupling because of how its ratings were
#' distributed, not because of how its factors were connected.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#' @param B Ensemble size. Draws that fail strong connectivity or hit the
#'   degenerate case are rejected and redrawn, so `B` rows are always returned.
#' @param seed Optional integer. Sets the random seed before drawing, so a
#'   figure can be reproduced. Passing `NULL` leaves the caller's random state
#'   alone.
#'
#' @return A data frame of `B` rows with columns `mu_max`, `lambda_max`,
#'   `dominance`, `hierarchy_sd`, `hierarchy_pr` and `hierarchy_gini`.
#'
#' @section Cost:
#' This is the one part of the package that is not instant. Each draw is an
#' eigendecomposition plus a connectivity check, so `B = 200` is two hundred of
#' each. At DEMATEL sizes that is still fast, but it is the only function here
#' that should be run on demand rather than eagerly.
#'
#' @section Rejection can loop:
#' A matrix sparse enough that most shuffles disconnect the graph will reject
#' many draws before accepting `B`. The loop has no iteration cap; this is
#' inherited from the source implementation and is noted in DEVELOPMENT.md.
#'
#' @examples
#' set.seed(1)
#' A <- matrix(sample(0:4, 36, replace = TRUE), 6, 6)
#' diag(A) <- 0
#' storage.mode(A) <- "double"
#' if (is_irreducible(A)) {
#'   ens <- surrogate_ensemble(A, B = 20, seed = 42)
#'   observed <- spectral_diagnostics(A)$mu_max
#'   mean(ens$mu_max >= observed)  # where the observation sits in the null
#' }
#'
#' @export
surrogate_ensemble <- function(A, B = 200, seed = NULL) {
  if (!is.null(seed)) set.seed(seed)
  n <- nrow(A)
  off <- which(row(A) != col(A))
  out <- vector("list", B)
  b <- 0
  while (b < B) {
    S <- A
    S[off] <- sample(A[off])
    if (!is_irreducible(S)) next
    # checks = FALSE for the return shape rather than for speed: this path needs
    # NULL so a bad draw can be rejected with `next`, where checks = TRUE would
    # return an NA-filled row. The saving is real but small -- measured at
    # roughly 0.05 s across B = 200 at 13 factors.
    d <- spectral_diagnostics(S, type = "A", checks = FALSE)
    if (is.null(d)) next
    b <- b + 1
    out[[b]] <- data.frame(mu_max = d$mu_max, lambda_max = d$lambda_max,
                           dominance = d$dominance, hierarchy_sd = d$hierarchy_sd,
                           hierarchy_pr = d$hierarchy_pr,
                           hierarchy_gini = d$hierarchy_gini)
  }
  do.call(rbind, out)
}
