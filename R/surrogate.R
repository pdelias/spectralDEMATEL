# The permutation surrogate baselines.
#
# The full shuffle is ported from R/spectral_dematel.R in the JDS paper
# repository. The shuffle within rows is ported from permute_rows() in that
# repository's R/baseline/spectral_baseline.R. Both draw in the order the
# paper's code draws, so a seed gives the same matrices in both places.

#' Permutation surrogate baseline
#'
#' Shuffles the off-diagonal entries of `A` and recomputes the diagnostics.
#' Repeating this gives a null distribution against which an observed
#' diagnostic can be read. Two schemes are available, and they hold different
#' things fixed.
#'
#' The question a null answers is the one a single number cannot: how much of
#' the structure is something the analyst built, and how much follows from the
#' ratings on their own. A system whose coupling sits in the middle of its own
#' surrogate ensemble has coupling because of how its ratings were distributed,
#' and not because of how its factors were connected.
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
#' @param scheme What the shuffle holds fixed. `"full"`, the default, permutes
#'   all off-diagonal entries, keeping the number of factors, the density and
#'   the multiset of ratings. `"rows"` permutes the off-diagonal entries inside
#'   each row, which also keeps every factor's dispatch total and the ratings
#'   it gave.
#'
#' @return A data frame with columns `mu_max`, `lambda_max`, `dominance`,
#'   `hierarchy_sd`, `hierarchy_pr`, `hierarchy_gini` and `hierarchy_h`,
#'   carrying four attributes: `requested` (the `B` asked for), `attempts`
#'   (shuffles drawn), `complete` (whether `B` rows were reached) and `scheme`.
#'   Returns `NULL` if no draw was admissible, if `A` cannot be shuffled, or if
#'   `B` is not a positive number, matching [spectral_diagnostics()] on
#'   inadmissible input.
#'
#'   **The row count is not guaranteed to be `B`.** Read `complete` before
#'   reporting anything derived from the ensemble.
#'
#'   Of the four hierarchy columns, `hierarchy_h`, `hierarchy_sd` and
#'   `hierarchy_gini` are **high when influence is concentrated** and
#'   `hierarchy_pr` is **low**.
#'
#' @section The two schemes ask different questions:
#' The full shuffle destroys everything except the ratings themselves. An
#' observation outside it says that the arrangement of the ratings matters. The
#' shuffle within rows keeps how much each factor dispatches and destroys only
#' where it sends it. An observation outside that null says that the targets of
#' influence matter, beyond the fact that some factors dispatch more than
#' others. A diagnostic can sit outside the first and inside the second.
#'
#' A seeded `"full"` ensemble draws the same shuffles as version 0.6.0.
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
#' So the loop is bounded, under both schemes. Hitting the bound is not an
#' error and does not throw: the function returns the draws it did get and says
#' so in `complete`, on the principle that a short null distribution the caller
#' knows about is more useful than no answer and far more useful than a hang.
#' `50 * B` tolerates an acceptance rate down to about 2% before it bites.
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
#'
#'   # The same question with every dispatch total held fixed:
#'   rows <- surrogate_ensemble(A, B = 20, seed = 42, scheme = "rows")
#'   mean(rows$mu_max >= observed)
#' }
#'
#' @export
surrogate_ensemble <- function(A, B = 200, seed = NULL, max_attempts = 50L * B,
                               scheme = c("full", "rows")) {
  scheme <- match.arg(scheme)
  if (!can_shuffle(A) || !is_count(B) || !is_count(max_attempts)) return(NULL)
  if (!is.null(seed)) set.seed(seed)
  out <- vector("list", B)
  b <- 0
  attempts <- 0
  while (b < B && attempts < max_attempts) {
    attempts <- attempts + 1
    S <- shuffle_matrix(A, scheme)
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
                           hierarchy_gini = d$hierarchy_gini,
                           hierarchy_h = d$hierarchy_h)
  }
  if (b == 0) return(NULL)
  res <- do.call(rbind, out[seq_len(b)])
  attr(res, "requested") <- B
  attr(res, "attempts")  <- attempts
  attr(res, "complete")  <- b >= B
  attr(res, "scheme")    <- scheme
  res
}

#' One shuffled copy of a matrix. The diagonal stays where it is.
#'
#' "full" permutes all off-diagonal entries together. "rows" permutes the
#' off-diagonal entries of each row among themselves, row by row from the
#' first, as permute_rows() in the paper's code does. The order matters: it is
#' what makes a seed give the same draws here and there.
#' @noRd
shuffle_matrix <- function(A, scheme) {
  S <- A
  if (scheme == "full") {
    off <- which(row(A) != col(A))
    S[off] <- sample(A[off])
    return(S)
  }
  n <- nrow(A)
  for (i in seq_len(n)) {
    idx <- seq_len(n)[-i]
    x <- A[i, idx]
    S[i, idx] <- x[sample.int(length(x))]
  }
  S
}

#' Can this be shuffled at all: a square numeric matrix of finite, non-negative
#' entries with at least two factors.
#' @noRd
can_shuffle <- function(A) {
  is.matrix(A) && is.numeric(A) && nrow(A) == ncol(A) && nrow(A) >= 2L &&
    all(is.finite(A)) && all(A >= 0)
}

#' One finite number of at least one, as an ensemble size has to be.
#' @noRd
is_count <- function(x) {
  is.numeric(x) && length(x) == 1L && is.finite(x) && x >= 1
}
