# spectralDEMATEL 0.3.1

* A `skipped` check no longer carries a value or a list of factors. Two call
  sites got this wrong independently -- `connectivity_margin` reported the
  factors hanging on a single link even when the graph was already
  disconnected, and `nonnegative` reported a count of zero on a matrix whose
  entries had never been examined. A number or a list of names beside a verdict
  that was never reached reads as evidence for it. The invariant is now
  enforced in one place, so no future call site can reintroduce it. Found while
  wiring the checks into the Shiny application.

* Check reasons read correctly when the count is one: "1 factor is" rather
  than "1 factors are". The count is frequently 1, since a single blank row is
  the commonest way to fail a check.

# spectralDEMATEL 0.3.0

* A vignette, `vignette("diagnosing-a-system")`, working through a full
  diagnosis on the two published matrices the package ships. Every number in it
  is computed at build time, so the documentation cannot drift from the code.

  It doubles as an acceptance test. Three claims in the first draft did not
  survive being checked against the output and were rewritten: the prominence
  reversal is real but milder than stated on this matrix, and the surrogate
  ensemble says something more interesting than expected — the observed
  coupling sits above all 200 draws and mode dominance below all 200, while
  hierarchy falls inside the ensemble on the low side.

* No changes to any function or definition.

# spectralDEMATEL 0.2.0

The assumption checks become data the engine returns, and nothing in the package
raises a condition any more for something a user can cause by pasting an odd
matrix.

## New

* `assumption_checks()` returns one row per condition, with a stable
  identifier to branch on, a verdict, a plain-language reason, the number
  behind it, and the factors at fault. The table has a fixed set of rows in a
  fixed order whatever is passed in, so a batch job can stack a hundred of them
  and get a rectangular result.

* Four verdicts, not two. `pass` and `fail` are the assumption holding or not.
  `warn` is the case a binary verdict cannot express — in scope, but close
  enough to a boundary that the reading is worth less than it looks. `skipped`
  means a prerequisite failed and the check was never evaluated; it is **not** a
  pass and must not be rendered as one.

* `spectral_diagnostics()` gains `checks` and `engine_version` components. The
  numbers no longer travel without the checks that say whether they mean
  anything, and every diagnosis records the version that produced it.

* An inadmissible matrix now returns the full diagnosis shape with every
  numeric field `NA`, plus the checks explaining why. A batch job over a
  hundred matrices never changes shape and never fails.

## Changed

* `dematel()` no longer calls `stopifnot()` or `warning()`. It returns `NULL`
  invisibly for anything it cannot process — not a numeric matrix, not square,
  missing or infinite or negative entries, entirely zero, or the uniform-totals
  degenerate case. `assumption_checks()` is where the reason lives.

* `spectral_diagnostics()` gains a `checks` argument. `checks = FALSE` returns
  the numbers alone, or `NULL` when they cannot be computed, and is what
  `surrogate_ensemble()` uses so a bad draw can be rejected.

* `is_irreducible()` and the `strong_connectivity` check now share a single
  reachability implementation, so the verdict and the factors it names cannot
  disagree.

* `dematel()` additionally rejects an all-zero matrix, which previously divided
  by zero.

## Thresholds

Two checks are continuous rather than binary. Neither is fitted, both are
recommendations, and both were chosen against measurements:

* **Coupling** warns at 0.95. The corpus of 117 published systems has a 99th
  percentile of 0.961 and a maximum of 0.965, so this flags roughly the top 3%.

* **Eigenvalue conditioning** warns at 5. Across 300 random admissible matrices
  of 5 to 30 factors, `ev_condition` ran from 1.0 to 2.0; a symmetric hub gives
  exactly 1, an asymmetric hub 3.2, and a long directed chain 28 to 33. Chains
  and deep hierarchies are where a first-order derivative stops being
  informative, and 5 sits in the gap between the two regimes.

## Note for anyone moving from `DEMATEL.Sensitivity`

That application computes several of these quantities independently and does not
agree with this package. Its mode dominance takes the subdominant eigenvalue by
real part rather than by modulus, which understates the ratio — by a factor of
three on one of the two worked matrices. Its `condition_number` is
`lambda_max / lambda_min`, a different quantity from `ev_condition` and one that
can be negative. Prefer this package's values.

# spectralDEMATEL 0.1.0

First release. The five working functions extracted from the research code
behind the source paper, ported unchanged, with a test suite.

## Included

* `dematel()`, `is_irreducible()`, `spectral_diagnostics()`,
  `surrogate_ensemble()`, `sensitivity_matrix()`.
* `worked_matrices`, the two transcribed systems from the source paper's
  Section 5.
* Dependencies: base R and `stats`. Nothing else.

## Behaviour deliberately not changed

The bodies were a verbatim port, including two things that contradicted the
engine's own boundary rules. Both were pinned by tests marked `INHERITED:` and
both are fixed in 0.2.0.
