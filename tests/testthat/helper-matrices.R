# Fixtures shared across the test files.
#
# Everything here is either a matrix with a spectrum someone can work out on
# paper, or a generator of admissible random matrices for the identity tests.

#' A random matrix that is admissible: square, non-negative, zero diagonal,
#' strongly connected, and not the uniform-totals degenerate case.
#'
#' Redraws until admissible, so it can be used inside a loop without guards.
random_admissible <- function(n, scale = 0:4, p = 0.6) {
  repeat {
    A <- matrix(sample(scale, n * n, replace = TRUE), n, n) *
      (matrix(stats::runif(n * n), n, n) < p)
    diag(A) <- 0
    storage.mode(A) <- "double"
    if (!is_irreducible(A)) next
    if (is.null(suppressWarnings(dematel(A)))) next
    return(A)
  }
}

# A pure 3-cycle. D^3 is a multiple of the identity, so the spectrum of D sits
# on a circle and every eigenvalue is available in closed form.
cycle3 <- matrix(c(0, 2, 0,
                   0, 0, 3,
                   4, 0, 0), 3, 3, byrow = TRUE)
storage.mode(cycle3) <- "double"

# Admissible, but T is defective: mu = -0.2 has algebraic multiplicity 2 and
# geometric multiplicity 1. Exercises the path where eigen() returns a
# non-diagonalisable decomposition.
defective <- matrix(c(0, 1, 0,
                      0, 0, 1,
                      2, 3, 0), 3, 3, byrow = TRUE)
storage.mode(defective) <- "double"

# Rank one, A = outer(a, b). D has one non-zero eigenvalue, so T does too and
# mode dominance must be exactly zero.
rank1_a <- c(1, 2, 3)
rank1_b <- c(1, 1, 2)
rank1 <- outer(rank1_a, rank1_b)
storage.mode(rank1) <- "double"

# The smallest interesting system.
two_by_two <- matrix(c(0, 2,
                       1, 0), 2, 2, byrow = TRUE)
storage.mode(two_by_two) <- "double"

# Two 2-cycles that never meet: fails strong connectivity.
block_diagonal <- matrix(0, 4, 4)
block_diagonal[1, 2] <- block_diagonal[2, 1] <- 1
block_diagonal[3, 4] <- block_diagonal[4, 3] <- 1
storage.mode(block_diagonal) <- "double"
