# Working on spectralDEMATEL

How the package is put together and what has to stay true. Read this before
changing anything in `R/`.

## The layout

```
R/dematel.R            dematel(), is_irreducible()
R/diagnostics.R        spectral_diagnostics(), gini() [internal]
R/surrogate.R          surrogate_ensemble()
R/sensitivity.R        sensitivity_matrix()
R/data.R               documentation for the worked_matrices dataset
data-raw/              the script that builds that dataset; not shipped
tests/testthat/        six test files, described below
```

Five exported functions, one internal helper, one dataset. Base R and `stats`,
nothing else, and that budget is not negotiable: a zero-dependency engine is the
cheapest thing to compile to WebAssembly and the easiest to trust.

## The rule the package exists to keep

**It computes and returns data. It never draws, never reads a file, never
prints, and never knows Shiny exists.**

The four that get violated first: no plotting or colour, no
`message()`/`cat()`/`print()` on a success path, no `stop()` for anything a user
can cause by pasting an odd matrix, and no file reading. Parsing belongs to
whatever is calling this.

Two of those are currently violated by inherited code — see *What was
deliberately left alone* below.

## The six test files and what each is for

They are not six ways of doing the same thing. Each catches a different class
of mistake, and dropping one leaves a real gap.

| File | Catches |
|---|---|
| `test-golden.R` | a definition changing at all |
| `test-known-answer.R` | a definition being *wrong* in a plausible-looking way |
| `test-identities.R` | two quantities drifting apart from each other |
| `test-sensitivity.R` | a closed form nobody differentiated numerically |
| `test-degenerate.R` | inputs the mathematics cannot take |
| `test-surrogate.R` | a baseline quietly becoming a different distribution |

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

## What was deliberately left alone

The bodies of the five functions are a **verbatim** port. Two things in them
contradict the rule above:

1. `dematel()` calls `stopifnot()` for non-square, negative and non-matrix
   inputs — all of which a user can produce by pasting.
2. The uniform-totals case signals through `warning()` plus a `NULL` return.
   A warning can only be caught and re-parsed, which is how message text becomes
   an accidental API.

Both are pinned by tests marked `INHERITED:` in `test-degenerate.R`. They are
step 3's job. Porting unchanged first means that when the assumption checks
land, the diff shows exactly what changed and the baseline it changed from was
tested.

One more inherited rough edge, not yet pinned: `surrogate_ensemble()` redraws
until it has `B` admissible shuffles, with no iteration cap. A matrix sparse
enough that most shuffles disconnect the graph will spin.

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
