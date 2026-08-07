# spectralDEMATEL 0.1.0

First release. The five working functions extracted from the JDS paper
repository's `R/spectral_dematel.R`, ported unchanged, with a test suite.

## Included

* `dematel()`, `is_irreducible()`, `spectral_diagnostics()`,
  `surrogate_ensemble()`, `sensitivity_matrix()`.
* `worked_matrices`, the two transcribed systems from the source paper's
  Section 5.
* Dependencies: base R and `stats`. Nothing else.

## Behaviour deliberately not changed

The bodies are a verbatim port. Two things in them contradict the engine's own
boundary rules and are scheduled for the next release rather than fixed here, so
that the change is a legible diff against a tested baseline:

* `dematel()` calls `stopifnot()` for conditions a user can cause by pasting an
  odd matrix (not square, negative entries, not a matrix).
* The uniform-totals case signals through `warning()` and a `NULL` return rather
  than a returned verdict.

Both are pinned by tests in `test-degenerate.R` under the heading `INHERITED:`.

## Note for anyone moving from `DEMATEL.Sensitivity`

That application computes several of these quantities independently and does not
agree with this package. Its mode dominance takes the subdominant eigenvalue by
real part rather than by modulus, which understates the ratio — by a factor of
three on one of the two worked matrices. Its `condition_number` is
`lambda_max / lambda_min`, a different quantity from `ev_condition` and one that
can be negative. Prefer this package's values.
