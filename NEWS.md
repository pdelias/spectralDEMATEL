# spectralDEMATEL 0.7.0

Not yet released. The development version is 0.6.0.9000, and this entry covers
the first milestone only.

Hierarchy is now classified on H, the reading of the revised source paper. This
changes the default of five functions. Read the first section before comparing
a result with one from 0.6.0.

## The default hierarchy reading changed

* `structural_type()`, `type_stability()`, `tradeoff_residual()`,
  `surrogate_position()` and `measurement_stability()` gain
  `reading = c("H", "sd")`. The default is `"H"`: the type is read on
  `hierarchy_h` with a cut of 0.10, against a reference corpus of 81 models.
  Up to 0.6.0 the type was read on `hierarchy_sd` with a cut of 0.10, against
  117 systems.

* **A matrix can receive another type than 0.6.0 gave it.** For the unit-length
  entry profile, H = (n - 1) sd(u)^2. A cut of 0.10 on sd is therefore a cut of
  0.10 on H only for eleven factors. A larger model can be hierarchical on H
  and diffuse on sd, and a smaller one the reverse. Both readings run the same
  way: high means concentrated.

* To reproduce a type from 0.6.0, name the old reading:

  ```r
  d <- spectral_diagnostics(A)
  structural_type(d, reading = "sd")
  ```

  The same argument reproduces `type_stability()`, `tradeoff_residual()`,
  `surrogate_position()` and `measurement_stability()`. Every value those
  functions returned in 0.6.0 is returned unchanged. The lists are longer, as
  described below.

* Every classification states the reading it used. `hierarchy_reading` is
  `"hierarchy_h"` or `"hierarchy_sd"` in the result of all five functions.
  `structural_type()` also returns `hierarchy_direction` and `corpus_n`.

## New

* `spectral_diagnostics()` returns `hierarchy_h`,
  H = 1 - (sum u)^2 / (n sum u^2), beside the three existing readings. It lies
  in [0, 1 - 1/n) and is **high when influence is concentrated**. A model whose
  factors are each split into copies keeps its H, while its `hierarchy_sd`
  falls. `surrogate_ensemble()` returns the same column.

* `structural_type()` returns `effort` under `reading = "H"`: the share of the
  attainable value that uniform effort forgoes (H), the share that appraisal on
  direct effects forgoes (mu_max^2), and the assumptions of the decision model
  in which both hold. The shares are `NA` under `reading = "sd"`.

* `surrogate_ensemble()` and `surrogate_position()` gain
  `scheme = c("full", "rows")`. `"rows"` shuffles the off-diagonal entries
  inside each row, which keeps every factor's dispatch total. A seeded
  `"full"` ensemble draws the same shuffles as 0.6.0.

* `surrogate_position()$metrics` gains a row for `hierarchy_h` and the columns
  `mean`, `sd` and `z`: the ensemble mean, its standard deviation, and the
  observed value less the mean in units of that deviation. `z` is `NA` when the
  ensemble does not vary.

* `normalising_line()` returns the normalising constant `s`, the rows and
  columns whose total equals it, and whether several lines tie.

* `prominence_relation()` returns dispatch, receipt, prominence, relation and
  the role of each factor, from the row and column totals of the total-relation
  matrix.

* `leverage_summary()` returns the line that sets `s`, the stabilisers, the
  amplifiers ranked by v_i u_j, and the eigenvalue condition number, in one
  object.

## Changed

* The internal `CORPUS` is now a list of two sets. `CORPUS$sd` holds the
  constants of 0.6.0, unchanged. `CORPUS$H` holds the aggregates of the
  81-model corpus: 57 diffuse-amplified, 13 hierarchical-dampened, 8
  hierarchical-amplified and 3 diffuse-dampened models, with median multipliers
  of 4.96, 1.54, 2.51 and 1.95. Its trade-off line is the median regression
  H = 0.374 - 0.436 mu_max, with a residual standard deviation of 0.075. Code
  that reads `spectralDEMATEL:::CORPUS` as a flat list has to choose a set.

* `sensitivity_matrix()` returns `line`, the result of `normalising_line()`, as
  a fourth component. `total`, `local` and `normalization` are unchanged.

* `type_stability()` takes `hierarchy_cuts = NULL` and then uses the grid of
  the reading: 0.05, 0.075, 0.10, 0.125 and 0.20 under `"H"`, and the grid of
  0.6.0 under `"sd"`.

* `tradeoff_residual()` returns `fit`, the method behind the corpus line:
  median regression under `"H"`, least squares under `"sd"`.

* `surrogate_ensemble()` carries a `scheme` attribute, and
  `surrogate_position()` returns `scheme`.

* Arguments that cannot be used return `NULL` and no longer raise an error: an
  ensemble size that is not a positive number in `surrogate_ensemble()`,
  `surrogate_position()` and `measurement_stability()`, and a cut that is not a
  finite number in `structural_type()` and `type_stability()`.
  `measurement_stability()` with `B = 0` returned a list of `NA` in 0.6.0 and
  returns `NULL` now.

## Where the package and the revised paper differ

* The paper zeroes the diagonal of a matrix before normalising. This package
  keeps the diagonal and reports it through the `zero_diagonal` check. A matrix
  with a non-zero diagonal therefore has a slightly different coupling here
  than in the paper. Every constant in `CORPUS$H` comes from matrices with a
  zero diagonal.

* For a model that is not strongly connected the paper admits shuffles that
  are not strongly connected either. `surrogate_position()` returns `NULL` for
  such a model, as in 0.6.0.

# spectralDEMATEL 0.6.0

The surrogate rejection loop is bounded.

## Changed

* `surrogate_ensemble()` gains `max_attempts`, default `50 * B`. It previously
  redrew until it had `B` admissible shuffles with no cap, so a matrix sparse
  enough that most shuffles disconnect the graph did not terminate in any useful
  time. Compiled to WebAssembly that is a frozen browser tab with no way to
  cancel, reachable by pasting a sparse matrix.

* **`surrogate_ensemble()` no longer guarantees `B` rows.** Reaching the bound
  returns the draws it did get rather than throwing, carrying `requested`,
  `attempts` and `complete` attributes. No admissible draw at all returns `NULL`,
  as `spectral_diagnostics()` does on inadmissible input.

* `surrogate_position()` gains `max_attempts` and returns `requested`,
  `complete` and `attempts` alongside `B`. **Read `complete` before quoting a
  surrogate share**: one computed over seven draws and one computed over two
  hundred are not the same claim, and `B` alone does not distinguish them.

A cap generous enough not to truncate reproduces the previous output exactly, so
ensembles drawn with an explicit seed against earlier versions are unchanged.

# spectralDEMATEL 0.5.0

The robustness section: whether a structural type is worth anything.

## New

* `surrogate_position()` places each observed diagnostic in its own permutation
  ensemble, holding the number of factors, the density and the exact multiset
  of ratings fixed. Reports the observed value, the ensemble range and median,
  the share of draws at or above the observation, and whether the observation
  falls outside the ensemble entirely. This is the comparison the source paper
  could not run across its corpus, because the precomputed spectral file holds
  no raw matrices.

* `measurement_stability()` perturbs the ratings by a stated tolerance and
  recomputes the type. Expert judgements on a four- or five-point scale carry
  at least half a point of noise, and a classification that does not survive
  that is not a classification. Recorded zeros are left alone by default: a
  zero is a judgement that there is no influence, not a one with noise on it.
  Perturbed values are clamped to the scale the analyst used.

Both take an exposed seed, so a figure that goes into a paper can be
reproduced.

## `type_share` is weak evidence for a common type

Nearly 70% of the reference corpus is diffuse-amplified, so a shuffled matrix
usually lands there too. A high `type_share` for a common type says the null
also produces that type, not that the observation is unremarkable. The
per-metric positions are the informative part. This is documented on the
function, because it is exactly the kind of number that reads as reassurance
when it carries none.

## A correction to how stability was tested

The first version of the test suite asserted that a larger tolerance never
leaves a type more secure than a smaller one. That was both vacuous and wrong.
Vacuous because it compared systems far from any boundary, where every
tolerance returns 1.000. Wrong because monotonicity does not hold in general: a
system sitting exactly on a cut stays near 0.5 at any tolerance, since more
noise pushes it across in both directions equally. It is replaced by a fixture
sitting one part in ten thousand from the hierarchy cut, which must report the
coin flip it is, and a moderate-distance case where monotonicity genuinely does
hold.

# spectralDEMATEL 0.4.0

The structure map: naming a system's type, and saying how firmly.

## New

* `structural_type()` places a diagnosed system on the (coupling, hierarchy)
  map, names its type, and returns the **margin to each cut** alongside. A
  system at coupling 0.51 and one at 0.95 are both amplified and the confidence
  in the advice is not remotely the same, so the margin travels with the name.
  It also carries the source paper's own wording for the intervention logic
  that type favours, and the caveat that the pairing is a hypothesis derived
  from structure rather than a validated result.

* `type_stability()` re-runs the classification across the range of hierarchy
  cuts anyone might reasonably choose, and reports whether the type survives.
  The cut is a recommendation, not a fitted constant, so this is the honest
  question rather than "which side is it on".

* `tradeoff_residual()` gives the system's distance from the corpus
  coupling-hierarchy trade-off, in residual standard deviations. A supplement
  to the computed hierarchy and never a substitute: used alone to assign types
  the trade-off is wrong for 19 of the 117 reference systems.

## Corpus constants

Aggregates only, in `CORPUS` (internal): the fitted trade-off with its residual
SD, counts and median multipliers by type, and 5/25/50/75/95 percentiles of
each axis. The per-system rows are deliberately **not** shipped -- they are
traceable to individual published studies and belong to work that is not yet
out. Everything the map and the residual reading need is here without them.

## The two cuts are not the same kind of constant

Coupling at 0.50 is the indirect-dominance threshold: at that value the
dominant eigenvalue of the total-relation matrix is exactly 1, so indirect
effects exactly equal direct ones. It follows from the algebra.

Hierarchy at 0.10 is a **recommendation**. `structural_type()` reports which
cuts it used and whether each was the default, so an interface can say so.

Only `hierarchy_sd` classifies. The 0.10 cut is calibrated against that reading
and every corpus-level result is expressed in it; `hierarchy_gini` and
`hierarchy_pr` measure the same idea on different scales, and the participation
ratio runs the opposite way. Show all three, classify on one.

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
