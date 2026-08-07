# Working on spectralDEMATEL

How the package is put together and what has to stay true. Read this before
changing anything in `R/`.

## The layout

```
R/dematel.R            dematel(), is_irreducible()
R/diagnostics.R        spectral_diagnostics(), diagnostics_core(), gini()
R/checks.R             assumption_checks(), reachability(), the thresholds
R/surrogate.R          surrogate_ensemble()
R/sensitivity.R        sensitivity_matrix()
R/data.R               documentation for the worked_matrices dataset
data-raw/              the script that builds that dataset; not shipped
tests/testthat/        seven test files, described below
```

Six exported functions, one dataset. Base R, `stats` and `utils`, nothing else,
and that budget is not negotiable: a near-zero-dependency engine is the cheapest
thing to compile to WebAssembly and the easiest to trust.

## The rule the package exists to keep

**It computes and returns data. It never draws, never reads a file, never
prints, and never knows Shiny exists.**

The four that get violated first: no plotting or colour, no
`message()`/`cat()`/`print()` on a success path, no `stop()` for anything a user
can cause by pasting an odd matrix, and no file reading. Parsing belongs to
whatever is calling this.

All four hold as of 0.2.0. `test-degenerate.R` is what keeps the third one true:
it runs every inadmissible input anyone could paste through every exported
function inside `expect_silent()`.

## The seven test files and what each is for

They are not seven ways of doing the same thing. Each catches a different class
of mistake, and dropping one leaves a real gap.

| File | Catches |
|---|---|
| `test-golden.R` | a definition changing at all |
| `test-known-answer.R` | a definition being *wrong* in a plausible-looking way |
| `test-identities.R` | two quantities drifting apart from each other |
| `test-sensitivity.R` | a closed form nobody differentiated numerically |
| `test-degenerate.R` | inputs the mathematics cannot take, and any condition raised for one |
| `test-surrogate.R` | a baseline quietly becoming a different distribution |
| `test-checks.R` | a ragged checks table, or `skipped` passing itself off as `pass` |

**The principle underneath all of them:** a metric definition has exactly one
implementation, and a test compares it against an *independent route* to the
same number. Agreement between a function and a literal it produced proves only
that nothing changed. The independent routes in use here are the corpus values
from the source paper's precomputed workbook (`test-golden.R`), closed-form
spectra worked out on paper (`test-known-answer.R`), the `type = "A"` and
`type = "T"` branches checked against each other (`test-identities.R`), and
finite differences (`test-sensitivity.R`).

**If a definition changes and the tests still pass, the tests were wrong.**

## Directions that invert, and have

Two things here read plausibly in either direction, and both have been wrong in
reviewed code in this project's history.

**Hierarchy has three readings running two ways.** `hierarchy_sd` and
`hierarchy_gini` are HIGH when influence enters at a few factors.
`hierarchy_pr` is LOW. The interface shows all three, each with a direction
label — that is a settled decision, not a default. `hierarchy_sd` is the reading
the source corpus uses throughout, so it is the one behind every corpus-level
comparison.

**Mode dominance takes the largest-modulus subdominant eigenvalue**, not the
largest real part. For a non-negative `T` the negative end of the spectrum
usually carries the larger modulus, so the real-part version understates the
ratio nearly everywhere while passing every plausibility check. It is wrong by a
factor of three on one of the two worked matrices, and it is what the currently
deployed application computes.

## The checks, and the two things that matter about them

`assumption_checks()` returns a fixed set of rows in a fixed order, whatever is
passed in. Two properties are load-bearing and both are tested:

**The table is rectangular.** A batch job stacks a hundred of them with `rbind`.
`CHECK_IDS` is both the declared order and the order the function emits, and
`skip_remaining()` indexes into it positionally — change one and you must change
the other. `test-checks.R` asserts they agree.

**`skipped` is not `pass`.** It means a prerequisite failed and the check was
never evaluated. An interface that renders the two the same way tells a user
their matrix is fine when nothing was tested. This is why there are four
verdicts rather than the three originally planned: a rectangular table needs a
row for a check that could not run.

The two continuous thresholds live in `CHECK_THRESHOLDS`, with the measurements
behind them in the comment above it. Neither is fitted. If either moves, say so
in `NEWS.md` — an interface quotes them to users as recommendations.

## Tolerances, and the one that bit twice

Most fixtures here have well-separated spectra and are asserted at 1e-10 to
1e-12, which holds on every platform. The `defective` fixture is different and
must not be treated the same way.

A defective eigenvalue — algebraic multiplicity 2, geometric multiplicity 1 —
perturbs as the *square root* of the backward error rather than linearly, so it
is computable to about `sqrt(.Machine$double.eps)`, roughly 1.5e-8, and no
further. macOS Accelerate and Linux LAPACK split the repeated pair differently.
An assertion at 1e-12 passed locally on macOS and failed on Linux; the corrected
1.49e-8 also failed, because `dominance` divides by `lambda_max` and amplifies
the absolute perturbation into a larger relative one.

The rule: **assert the perturbation bound tightly, and the quantity derived from
it loosely.** Anything reading a repeated or near-repeated eigenvalue gets a
tolerance derived from `sqrt(eps)` and the amplification, with the derivation
written down. A local macOS run cannot catch this — the Linux jobs in the CI
matrix are the only thing that can.

## Known rough edge

`surrogate_ensemble()` redraws until it has `B` admissible shuffles, with no
iteration cap. A matrix sparse enough that most shuffles disconnect the graph
will spin. Inherited from the source implementation, not yet pinned by a test.

## Changing something

1. **A metric definition** — ask first. Then change the one implementation,
   update the golden literals, bump the version, and write it in `NEWS.md`. Old
   diagnoses must stay interpretable rather than silently mean something new.
2. **A new function** — it returns data, and it gets tests from at least two of
   the routes above before it is exported.
3. **A dependency** — ask first. It must be in the webR repository, and the
   engine's answer is almost always no.

```bash
Rscript -e 'roxygen2::roxygenise(); devtools::test()'
```

Before anything is pushed:

```bash
R CMD build spectralDEMATEL && R CMD check --as-cran spectralDEMATEL_*.tar.gz
```

Currently `Status: OK` — no errors, no warnings, no notes. Keep it there; the
r-universe build that produces the WebAssembly binary is less forgiving than a
local `devtools::load_all()`, and the application cannot see a change to this
package until that binary exists.

## Regenerating the dataset

```bash
Rscript data-raw/worked_matrices.R
```

Reads `../reference/matrices/`. Both matrices carry `doi` and `label`
attributes. If either is ever re-transcribed, the golden literals in
`test-golden.R` must be re-derived from the corpus, not from this package.
