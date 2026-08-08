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
#' @param B Ensemble size requested. Draws that fail strong connectivity or hit
#'   the degenerate case are rejected and redrawn.
#' @param seed Optional integer. Sets the random seed before drawing, so a
#'   figure can be reproduced. Passing `NULL` leaves the caller's random state
#'   alone.
#' @param max_attempts Maximum number of shuffles to draw before giving up,
#'   default `50 * B`. Reaching it returns a short ensemble rather than
#'   continuing; see the section below.
#'
#' @return A data frame with columns `mu_max`, `lambda_max`, `dominance`,
#'   `hierarchy_sd`, `hierarchy_pr` and `hierarchy_gini`, carrying three
#'   attributes: `requested` (the `B` asked for), `attempts` (shuffles drawn)
#'   and `complete` (whether `B` rows were reached). Returns `NULL` if no draw
#'   was admissible, matching [spectral_diagnostics()] on inadmissible input.
#'
#'   **The row count is not guaranteed to be `B`.** Read `complete` before
#'   reporting anything derived from the ensemble.
#'
#' @section Cost:
#' This is the one part of the package that is not instant. Each draw is an
#' eigendecomposition plus a connectivity check, so `B = 200` is two hundred of
#' each. At DEMATEL sizes that is still fast, but it is the only function here
#' that should be run on demand rather than eagerly.
#'
#' @section When rejection dominates:
#' A matrix sparse enough that most shuffles disconnect the graph rejects most
#' draws. A directed cycle is the extreme case: its edges are exactly as many as
#' strong connectivity needs, so almost every rearrangement of them breaks it.
#'
#' The original implementation looped until it had `B` admissible draws with no
#' cap, which on such a matrix does not terminate in any useful time. In a
#' script that is a long wait; compiled to WebAssembly it is a frozen browser
#' tab with no way to cancel, reachable by pasting a sparse matrix.
#'
#' So the loop is bounded. Hitting the bound is not an error and does not throw:
#' the function returns the draws it did get and says so in `complete`, on the
#' principle that a short null distribution the caller knows about is more
#' useful than no answer and far more useful than a hang. `50 * B` tolerates an
#' acceptance rate down to about 2% before it bites.
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
surrogate_ensemble <- function(A, B = 200, seed = NULL, max_attempts = 50L * B) {
  if (!is.null(seed)) set.seed(seed)
  off <- which(row(A) != col(A))
  out <- vector("list", B)
  b <- 0
  attempts <- 0
  while (b < B && attempts < max_attempts) {
    attempts <- attempts + 1
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
  if (b == 0) return(NULL)
  res <- do.call(rbind, out[seq_len(b)])
  attr(res, "requested") <- B
  attr(res, "attempts")  <- attempts
  attr(res, "complete")  <- b >= B
  res
}
