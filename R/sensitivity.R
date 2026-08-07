# Closed-form derivative of the dominant eigenvalue in every off-diagonal entry.
#
# Ported verbatim from R/spectral_dematel.R in the JDS paper repository, where it
# is checked against finite differences by verify_spectral_claims.R. That check
# now lives in tests/testthat/test-sensitivity.R.

#' Per-link sensitivity of the dominant eigenvalue
#'
#' Returns \eqn{\partial \lambda_{max} / \partial a_{ij}} for every off-diagonal
#' entry of the direct influence matrix, split into the two effects that
#' produce it.
#'
#' This moves the question from "which factor matters" to "which relationship,
#' if strengthened or weakened, moves the system most", which is much closer to
#' what an intervention actually is. The derivative is signed: entries with a
#' negative derivative *lower* coupling when strengthened, and are the levers
#' for a system someone wants to calm rather than drive.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#'
#' @return A list of three `n` by `n` matrices, or `NULL` when [dematel()] found
#'   the degenerate case. Diagonals are `NA` throughout, since a self-loop is
#'   not an actionable link.
#'
#'   \describe{
#'     \item{`local`}{The effect of the entry on the propagation structure,
#'       holding the normalising constant fixed. Always positive.}
#'     \item{`normalization`}{The effect of the entry on the normalising
#'       constant `s`. Non-zero only in the row or column that currently sets
#'       `s`, and always negative there.}
#'     \item{`total`}{Their sum, which is the quantity to rank on.}
#'   }
#'
#' @section Why the split is not decoration:
#' `s` is the larger of the maximum row sum and the maximum column sum.
#' Increasing an entry that sits in the row setting `s` raises `s` too, which
#' shrinks the whole normalised matrix. The two effects oppose each other, and
#' in that row the normalisation effect can win — so raising an influence
#' rating there *lowers* the dominant eigenvalue. A user shown only the total
#' will read that as a bug. Showing the components names it as a mechanism.
#'
#' @section Never show this without the condition number:
#' These are first-order estimates. `ev_condition` from
#' [spectral_diagnostics()] bounds how far they can be trusted, and large values
#' mean the derivative is locally uninformative. The ranking without it
#' misleads, so the two ship together or neither ships.
#'
#' @examples
#' A <- matrix(c(0, 2, 1,
#'               3, 0, 1,
#'               1, 2, 0), 3, 3, byrow = TRUE)
#' s <- sensitivity_matrix(A)
#'
#' # The single most effective link to strengthen:
#' which(s$total == max(s$total, na.rm = TRUE), arr.ind = TRUE)
#'
#' # How far to trust that: 1 is ideal, large is uninformative.
#' spectral_diagnostics(A)$ev_condition
#'
#' @export
sensitivity_matrix <- function(A) {
  m <- dematel(A); if (is.null(m)) return(NULL)
  n <- nrow(A); s <- m$s; T <- m$T
  eT <- eigen(T);  k  <- which.max(Re(eT$values));  u <- Re(eT$vectors[, k])
  eTt <- eigen(t(T)); kl <- which.max(Re(eTt$values)); v <- Re(eTt$vectors[, kl])
  v <- v / sum(v * u)
  S <- diag(n) + T
  C <- as.numeric(t(v) %*% S %*% A %*% S %*% u)
  local <- outer(as.vector(t(v) %*% S), as.vector(S %*% u)) / s
  rs <- rowSums(A); cs <- colSums(A)
  norm_eff <- matrix(0, n, n)
  if (max(rs) >= max(cs)) norm_eff[abs(rs - s) < 1e-12, ] <- -C / s^2
  if (max(cs) >= max(rs)) norm_eff[, abs(cs - s) < 1e-12] <- -C / s^2
  total <- local + norm_eff
  diag(local) <- diag(norm_eff) <- diag(total) <- NA  # diagonal not actionable
  list(total = total, local = local, normalization = norm_eff)
}
