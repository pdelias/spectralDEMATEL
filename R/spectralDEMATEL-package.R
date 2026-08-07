#' spectralDEMATEL: spectral diagnostics for DEMATEL influence matrices
#'
#' The premise of the method these functions implement is that you cannot look
#' up your intervention logic, you have to diagnose it. Two DEMATEL models with
#' the same number of factors and the same rating scale can require opposite
#' interventions, and which one you are holding is a property of the matrix's
#' spectrum rather than of its subject matter.
#'
#' @section The five working functions:
#' \describe{
#'   \item{[dematel()]}{Normalisation and the total-relation matrix.}
#'   \item{[is_irreducible()]}{Strong connectivity, assumption A2.}
#'   \item{[spectral_diagnostics()]}{The reported diagnostic set.}
#'   \item{[surrogate_ensemble()]}{A null distribution preserving size, density
#'     and the entry distribution.}
#'   \item{[sensitivity_matrix()]}{Per-link derivative of the dominant
#'     eigenvalue.}
#' }
#'
#' @section What this package does not do:
#' It computes and returns data. It does not draw, does not read files, does not
#' print on a success path, and does not know Shiny exists. Parsing a
#' spreadsheet, choosing a colour and rendering a verdict all belong to whatever
#' is calling it. This is what lets the same code serve a script, a batch job
#' over a hundred matrices, and a browser-side application compiled to
#' WebAssembly, without a second implementation.
#'
#' @section Reading the numbers in the right direction:
#' Two definitions here invert easily and have been inverted before, in reviewed
#' code, while producing entirely plausible output.
#'
#' Hierarchy has three readings and they do not run the same way.
#' `hierarchy_sd` and `hierarchy_gini` are **high** when influence enters at a
#' few factors; `hierarchy_pr` is **low**. Anything displaying one of them must
#' say which reading and which direction.
#'
#' Mode dominance takes the **largest-modulus** subdominant eigenvalue, not the
#' largest real part. For a non-negative total-relation matrix the negative end
#' of the spectrum usually carries the larger modulus.
#'
#' @section One quantity, three scales:
#' `mu_max` is coupling on \eqn{(0, 1)}; `lambda_max` is
#' \eqn{\mu/(1 - \mu)}; and the multiplier \eqn{1/(1 - \mu_{max})} equals
#' \eqn{1 + \lambda_{max}}. They are the same fact reported three ways, and the
#' tests assert the identities rather than assuming them.
#'
#' @keywords internal
"_PACKAGE"
