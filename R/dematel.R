# The DEMATEL pipeline: normalisation, the total-relation matrix, and the
# strong-connectivity check.
#
# Ported verbatim from R/spectral_dematel.R in the JDS paper repository, which
# is the single source of truth for every metric definition here. The bodies
# below are unchanged from that file. See ../reference/engine/ for the original
# and DEVELOPMENT.md for what "verbatim" is protecting.

#' Normalise a direct influence matrix and compute the total-relation matrix
#'
#' Divides `A` by the larger of its maximum row sum and its maximum column sum
#' to obtain the normalised matrix `D`, then sums the influence propagated over
#' paths of every length to obtain the total-relation matrix
#' \eqn{T = D(I - D)^{-1}}.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#'
#' @return A list with components `A` (the input), `s` (the normalising
#'   constant), `D` (the normalised matrix) and `T` (the total-relation
#'   matrix). Returns `NULL` invisibly, with a warning, when the input is
#'   degenerate.
#'
#' @section The degenerate case:
#' When every row and column total is the same, the spectral radius of `D` is
#' exactly 1, \eqn{I - D} is singular and `T` is undefined (Lee et al. 2013).
#' This is a property of the input rather than a mistake by the caller, and a
#' uniform matrix is the commonest way to meet it. The function warns and
#' returns `NULL`; every function here that calls it propagates the `NULL`
#' rather than failing.
#'
#' @section Not yet an assumption check:
#' The `stopifnot()` below and the `warning()` above are inherited from the
#' paper's implementation and are scheduled for replacement by returned
#' verdicts. Do not build on them. See DEVELOPMENT.md.
#'
#' @examples
#' A <- matrix(c(0, 2, 0,
#'               0, 0, 3,
#'               4, 0, 0), 3, 3, byrow = TRUE)
#' m <- dematel(A)
#' m$s
#' round(m$T, 4)
#'
#' @export
dematel <- function(A) {
  stopifnot(is.matrix(A), nrow(A) == ncol(A), all(A >= 0))
  s <- max(max(rowSums(A)), max(colSums(A)))
  D <- A / s
  mu_max <- max(Re(eigen(D, only.values = TRUE)$values))
  if (mu_max >= 1 - 1e-12) {
    warning("rho(D) = 1: uniform-totals pathology, T undefined (Lee et al. 2013).")
    return(invisible(NULL))
  }
  T <- D %*% solve(diag(nrow(A)) - D)
  list(A = A, s = s, D = D, T = T)
}

#' Strong connectivity of the influence graph (assumption A2)
#'
#' Tests whether the directed graph whose edges are the positive entries of `A`
#' is strongly connected, that is, whether every factor can reach every other
#' factor along some directed path. This is assumption A2, and it is what makes
#' the dominant eigenvector unique and strictly positive.
#'
#' The test is performed by repeatedly squaring the boolean reachability matrix
#' \eqn{(I + M)}, where `M` is the sign pattern of `A`. Costs \eqn{O(n^4)} in
#' the worst case, which is free at the sizes DEMATEL models reach.
#'
#' @param A A square numeric matrix. Only the sign pattern is used.
#'
#' @return `TRUE` if the influence graph is strongly connected, `FALSE`
#'   otherwise.
#'
#' @section Failure is common and is not fatal:
#' A single factor that dispatches nothing gives an all-zero row and fails the
#' test. One of the two worked matrices in the source paper fails it, and its
#' diagnostics are reported regardless. Callers should report the verdict rather
#' than refuse to compute.
#'
#' @examples
#' cycle <- matrix(c(0, 1, 0,
#'                   0, 0, 1,
#'                   1, 0, 0), 3, 3, byrow = TRUE)
#' is_irreducible(cycle)
#'
#' # A factor that influences nothing breaks strong connectivity.
#' sink_node <- cycle
#' sink_node[3, ] <- 0
#' is_irreducible(sink_node)
#'
#' @export
is_irreducible <- function(A) {
  n <- nrow(A)
  M <- (A > 0) * 1
  R <- diag(n) + M
  P <- R
  for (k in seq_len(n)) P <- (P %*% R > 0) * 1
  all(P > 0)
}
