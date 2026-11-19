# The line that sets the normalising constant, prominence and relation, and the
# leverage summary.
#
# Three views of where a system can be moved. None of them introduces a new
# definition: the line is the rule sensitivity_matrix() has always applied, the
# totals are the row and column sums of the total-relation matrix, and the
# summary arranges what sensitivity_matrix() and spectral_diagnostics() return.

#' The row or column that sets the normalising constant
#'
#' DEMATEL divides the direct influence matrix by `s`, the largest of its row
#' and column totals. This function reports `s` and which rows and columns
#' attain it. Those lines matter for intervention: a rating on one of them
#' raises `s` when it is raised, which shrinks the whole normalised matrix, and
#' it is only there that strengthening a link can lower the coupling of the
#' system.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#'
#' @return A list, or `NULL` when `A` is not a square numeric matrix of finite,
#'   non-negative entries with a positive total:
#'   \describe{
#'     \item{`s`}{The normalising constant.}
#'     \item{`rows`}{Indices of the rows whose total equals `s`.}
#'     \item{`cols`}{Indices of the columns whose total equals `s`.}
#'     \item{`tied`}{`TRUE` when more than one line attains `s`.}
#'   }
#'
#' @section When lines tie:
#' With one line the derivative of the dominant eigenvalue in a rating on that
#' line is well defined. With several, raising a rating on one of them raises
#' `s` and lowering it does not, because another line still holds the maximum.
#' The sensitivities on a tied line therefore describe an increase only. A tie
#' is common in small matrices of integer ratings, and `tied` is returned so
#' that it can be reported beside the ranking.
#'
#' A total counts as equal to `s` when it lies within `1e-12` of it, the rule
#' [sensitivity_matrix()] has applied since its first version.
#'
#' @section It does not need the total-relation matrix:
#' When every row and column total is the same, [dematel()] returns `NULL`
#' because the total-relation matrix is undefined. The line is still a fact
#' about such a matrix, and this function returns it, with every line tied.
#'
#' @examples
#' A <- matrix(c(0, 3, 3,
#'               1, 0, 1,
#'               1, 1, 0), 3, 3, byrow = TRUE)
#' normalising_line(A)
#'
#' @seealso [sensitivity_matrix()], which returns the same object as its
#'   `line` component, and [leverage_summary()].
#' @export
normalising_line <- function(A) {
  if (!is.matrix(A) || !is.numeric(A)) return(NULL)
  if (nrow(A) != ncol(A) || nrow(A) < 1L) return(NULL)
  if (!all(is.finite(A)) || any(A < 0)) return(NULL)

  rs <- unname(rowSums(A))
  cs <- unname(colSums(A))
  s <- max(max(rs), max(cs))
  if (!is.finite(s) || s <= 0) return(NULL)

  # Ported from sensitivity_matrix(), where the same two conditions chose the
  # entries that receive the normalisation effect.
  rows <- if (max(rs) >= max(cs)) which(abs(rs - s) < 1e-12) else integer(0)
  cols <- if (max(cs) >= max(rs)) which(abs(cs - s) < 1e-12) else integer(0)

  list(s = s, rows = rows, cols = cols,
       tied = length(rows) + length(cols) > 1L)
}

#' Prominence and relation of each factor
#'
#' The standard product of a DEMATEL study: the row and column totals of the
#' total-relation matrix, their sum and their difference.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#'
#' @return A data frame with one row per factor, or `NULL` when [dematel()]
#'   cannot process `A`:
#'   \describe{
#'     \item{`factor`}{The row name of `A`, or the column name, or the index of
#'       the factor as text when `A` does not have names.}
#'     \item{`dispatch`}{Row total of the total-relation matrix: the influence
#'       the factor sends, directly and indirectly.}
#'     \item{`receipt`}{Column total: the influence it receives.}
#'     \item{`prominence`}{`dispatch + receipt`.}
#'     \item{`relation`}{`dispatch - receipt`.}
#'     \item{`role`}{`"net cause"` when `relation` is positive, `"net effect"`
#'       when it is negative, and `NA` when it is exactly zero.}
#'   }
#'
#' @section What these do not say:
#' Prominence adds what a factor sends to what it receives, so a factor that
#' gathers influence can outrank the one that injects it. The entry and
#' accumulation profiles of [spectral_diagnostics()] keep the two roles apart.
#' A relation close to zero gives a role that a small change in one rating can
#' reverse, so the sign is best read with the size.
#'
#' @examples
#' A <- matrix(c(0, 2,
#'               1, 0), 2, 2, byrow = TRUE)
#' prominence_relation(A)
#'
#' @seealso [spectral_diagnostics()] for the entry and accumulation profiles.
#' @export
prominence_relation <- function(A) {
  m <- dematel(A)
  if (is.null(m)) return(NULL)

  dispatch <- unname(rowSums(m$T))
  receipt  <- unname(colSums(m$T))
  relation <- dispatch - receipt

  role <- rep(NA_character_, length(relation))
  role[relation > 0] <- "net cause"
  role[relation < 0] <- "net effect"

  data.frame(
    factor     = factor_labels(A),
    dispatch   = dispatch,
    receipt    = receipt,
    prominence = dispatch + receipt,
    relation   = relation,
    role       = role,
    stringsAsFactors = FALSE)
}

#' Labels for the factors of a matrix: row names, else column names, else the
#' index as text.
#' @noRd
factor_labels <- function(A) {
  if (!is.null(rownames(A))) return(rownames(A))
  if (!is.null(colnames(A))) return(colnames(A))
  as.character(seq_len(nrow(A)))
}

#' Where the system can be moved, with what bounds the claim
#'
#' Arranges the per-link sensitivities of [sensitivity_matrix()] into the two
#' groups an intervention reads them in, and returns the eigenvalue condition
#' number in the same object, so that the ranking cannot be reported without
#' it.
#'
#' @param A A square, non-negative numeric matrix of direct influences.
#'
#' @return A list, or `NULL` when [dematel()] cannot process `A`:
#'   \describe{
#'     \item{`line`}{The result of [normalising_line()].}
#'     \item{`stabilisers`}{Data frame of the links whose total derivative is
#'       negative, the most negative first: `from`, `to` (indices of the
#'       dispatching and the receiving factor), `rating` (the current entry of
#'       `A`), `total`, `local` and `normalization` as in
#'       [sensitivity_matrix()]. Strengthening such a link lowers the dominant
#'       eigenvalue. Every one lies on `line`.}
#'     \item{`amplifiers`}{Data frame of the links whose total derivative is
#'       positive: `from`, `to`, `rating`, `weight` and `total`, ranked by
#'       `weight`, the largest first. `weight` is \eqn{v_i u_j}, the
#'       accumulation weight of the dispatching factor times the entry weight
#'       of the receiving one, with the two profiles scaled so that
#'       \eqn{v'u = 1}.}
#'     \item{`ev_condition`}{The eigenvalue condition number of
#'       [spectral_diagnostics()]. At least 1. Large values mean the
#'       derivatives are locally uninformative.}
#'     \item{`mu_max`}{Coupling. A link on the line is a stabiliser exactly
#'       when its `weight` is below it.}
#'     \item{`engine_version`}{The package version that produced the summary.}
#'   }
#'
#' @section Reading the two tables:
#' Off the line that sets `s`, every link is an amplifier, and its derivative is
#' proportional to `weight`. Ranking by `weight` therefore ranks the links off
#' the line by their effect. On the line a second term enters with the opposite
#' sign, and the link is a stabiliser when that term wins. The rule is a
#' comparison of two numbers: a link on the line is a stabiliser when its
#' `weight` is below `mu_max`.
#'
#' A link with a `rating` of zero is listed like any other. The derivative says
#' what would follow if an influence were introduced there.
#'
#' @section First-order, and bounded:
#' These are derivatives, so they describe a small change. `ev_condition`
#' bounds how far they can be trusted. When `line$tied` is `TRUE` the
#' derivative on the tied lines describes an increase of the rating only.
#'
#' @examples
#' A <- matrix(c(0, 2, 1,
#'               3, 0, 1,
#'               1, 2, 0), 3, 3, byrow = TRUE)
#' l <- leverage_summary(A)
#' l$line
#' l$stabilisers
#' head(l$amplifiers, 3)
#' l$ev_condition
#'
#' @seealso [sensitivity_matrix()] for the full matrices.
#' @export
leverage_summary <- function(A) {
  sens <- sensitivity_matrix(A)
  if (is.null(sens)) return(NULL)
  d <- spectral_diagnostics(A, checks = FALSE)
  if (is.null(d)) return(NULL)

  # v_i u_j with v'u = 1, from the two profiles of the diagnosis.
  weight <- outer(d$accumulation, d$entry_points) /
    sum(d$accumulation * d$entry_points)

  links <- which(row(A) != col(A), arr.ind = TRUE)
  at <- function(M) unname(M[links])
  all_links <- data.frame(
    from          = unname(links[, 1]),
    to            = unname(links[, 2]),
    rating        = at(A),
    weight        = at(weight),
    total         = at(sens$total),
    local         = at(sens$local),
    normalization = at(sens$normalization),
    stringsAsFactors = FALSE)

  stab <- all_links[all_links$total < 0, , drop = FALSE]
  stab <- stab[order(stab$total, stab$from, stab$to),
               c("from", "to", "rating", "total", "local", "normalization"),
               drop = FALSE]
  rownames(stab) <- NULL

  amp <- all_links[all_links$total > 0, , drop = FALSE]
  amp <- amp[order(-amp$weight, amp$from, amp$to),
             c("from", "to", "rating", "weight", "total"), drop = FALSE]
  rownames(amp) <- NULL

  list(
    line           = sens$line,
    stabilisers    = stab,
    amplifiers     = amp,
    ev_condition   = d$ev_condition,
    mu_max         = d$mu_max,
    engine_version = engine_version()
  )
}
