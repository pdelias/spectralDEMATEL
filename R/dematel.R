# The DEMATEL pipeline: normalisation, the total-relation matrix, and the
# strong-connectivity check.
#
# The arithmetic is the paper's R/spectral_dematel.R unchanged. What changed in
# 0.2.0 is the failure path: these functions no longer stop() or warning() for
# anything a user can cause by pasting an odd matrix. They return NULL, and
# assumption_checks() explains why.

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
#'   matrix). Returns `NULL` invisibly for any input it cannot process.
#'
#' @section Why it returns NULL instead of failing:
#' A user pasting an asymmetric matrix is not an exceptional condition, it is
#' Tuesday. Every rejection here is silent and recoverable, and
#' [assumption_checks()] is where the reason lives — as data, in plain language,
#' naming the factors at fault. Nothing in this package raises a condition for
#' something a user can produce.
#'
#' `NULL` is returned when `A` is not a numeric matrix, is not square, holds a
#' missing, infinite or negative entry, is entirely zero, or meets the
#' degenerate case below.
#'
#' @section The degenerate case:
#' When every row and column total is the same, the spectral radius of `D` is
#' exactly 1, \eqn{I - D} is singular and `T` is undefined (Lee et al. 2013).
#' A matrix in which every pair was rated identically always meets it.
#'
#' @examples
#' A <- matrix(c(0, 2, 0,
#'               0, 0, 3,
#'               4, 0, 0), 3, 3, byrow = TRUE)
#' m <- dematel(A)
#' m$s
#' round(m$T, 4)
#'
#' # Degenerate input returns NULL rather than failing:
#' is.null(dematel(matrix(2, 4, 4)))
#'
#' @seealso [assumption_checks()], which reports why an input was rejected.
#' @export
dematel <- function(A) {
  if (!is.matrix(A) || !is.numeric(A)) return(invisible(NULL))
  if (nrow(A) != ncol(A) || nrow(A) < 1L) return(invisible(NULL))
  if (!all(is.finite(A)) || any(A < 0)) return(invisible(NULL))

  s <- max(max(rowSums(A)), max(colSums(A)))
  if (!is.finite(s) || s <= 0) return(invisible(NULL))

  D <- A / s
  mu_max <- max(Re(eigen(D, only.values = TRUE)$values))
  if (!is.finite(mu_max) || mu_max >= 1 - 1e-12) return(invisible(NULL))

  T <- tryCatch(D %*% solve(diag(nrow(A)) - D), error = function(e) NULL)
  if (is.null(T) || !all(is.finite(T))) return(invisible(NULL))

  list(A = A, s = s, D = D, T = T)
}

#' Strong connectivity of the influence graph (assumption A2)
#'
#' Tests whether the directed graph whose edges are the positive entries of `A`
#' is strongly connected, that is, whether every factor can reach every other
#' factor along some directed path. This is assumption A2, and it is what makes
#' the dominant eigenvector unique and strictly positive.
#'
#' @param A A square numeric matrix. Only the sign pattern is used.
#'
#' @return `TRUE` if the influence graph is strongly connected, `FALSE`
#'   otherwise.
#'
#' @section Failure is common and is not fatal:
#' A single factor that dispatches nothing gives an all-zero row and fails the
#' test. One of the two matrices in [worked_matrices] fails it, and its
#' diagnostics are published regardless. Use
#' `assumption_checks()` to learn *which* factors are stranded.
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
  all(reachability(A) > 0)
}
