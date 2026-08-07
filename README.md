# spectralDEMATEL

<!-- badges: start -->
[![R-CMD-check](https://github.com/pdelias/spectralDEMATEL/actions/workflows/R-CMD-check.yaml/badge.svg)](https://github.com/pdelias/spectralDEMATEL/actions/workflows/R-CMD-check.yaml)
[![r-universe](https://pdelias.r-universe.dev/badges/spectralDEMATEL)](https://pdelias.r-universe.dev/spectralDEMATEL)
<!-- badges: end -->

Spectral diagnostics for DEMATEL influence matrices.

Two DEMATEL models with the same number of factors and the same rating scale can
require opposite interventions. Which one you are holding is a property of the
matrix's spectrum, not of its subject matter — so **you cannot look up your
intervention logic, you have to diagnose it.** This package is the diagnosis.

## Install

```r
install.packages("spectralDEMATEL",
                 repos = c("https://pdelias.r-universe.dev",
                           "https://cloud.r-project.org"))
```

In [webR](https://docs.r-wasm.org/webr/latest/) or a
[shinylive](https://posit-dev.github.io/r-shinylive/) application, point at the
WebAssembly builds instead:

```r
install.packages("spectralDEMATEL",
                 repos = c("https://pdelias.r-universe.dev",
                           "https://repo.r-wasm.org"))
```

## Use

```r
library(spectralDEMATEL)

A <- worked_matrices$fefo_stock_control   # 13 factors, from a published study
d <- spectral_diagnostics(A)

d$mu_max              # 0.828  coupling, in (0, 1)
d$multiplier          # 5.81   total effects are 5.8x the direct ones
d$indirect_dominant   # TRUE   indirect effects exceed direct ones
d$dominance           # 0.028  one propagation mode, no competing second
d$hierarchy_sd        # 0.038  low, so influence enters diffusely
d$ev_condition        # 1.01   sensitivity estimates are trustworthy here
```

Where to intervene:

```r
s <- sensitivity_matrix(A)
which(s$total == max(s$total, na.rm = TRUE), arr.ind = TRUE)
```

`s$total` is the derivative of the dominant eigenvalue in every off-diagonal
entry, so it ranks *relationships* rather than factors. It is signed: negative
entries are the levers for calming a system rather than driving it. Read it
alongside `d$ev_condition`, which bounds how far a first-order estimate can be
trusted — **never one without the other**.

Is the structure real, or does it follow from the entry distribution alone?

```r
ens <- surrogate_ensemble(A, B = 200, seed = 42)
mean(ens$mu_max >= d$mu_max)
```

## Reading the numbers in the right direction

Two definitions here invert easily, and both have been inverted before in
reviewed code, while producing entirely plausible output.

**Hierarchy has three readings and they do not run the same way.**

| | High means | Size-free |
|---|---|---|
| `hierarchy_sd` | concentrated | no |
| `hierarchy_gini` | concentrated | yes |
| `hierarchy_pr` | **diffuse** — the opposite | yes |

Anything that displays one of these must say which reading it used and which
way it runs. `hierarchy_sd` is the reading behind every published corpus result.

**Mode dominance takes the largest-modulus subdominant eigenvalue**, not the
largest real part. For a non-negative total-relation matrix the negative end of
the spectrum usually carries the larger modulus, so the real-part version
understates the ratio nearly everywhere — while passing every plausibility
check.

**Coupling, the dominant eigenvalue and the multiplier are one quantity on three
scales:** `mu_max` in (0, 1), `lambda_max = mu/(1 - mu)`, and
`multiplier = 1/(1 - mu_max) = 1 + lambda_max`.

## What it does not do

It computes and returns data. It does not draw, read files, print on a success
path, or know Shiny exists. Parsing a spreadsheet, choosing a colour and
rendering a verdict belong to whatever is calling it.

That boundary is what lets the same code serve a script, a batch job over a
hundred matrices, and a browser-side application compiled to WebAssembly,
without a second implementation — and a second implementation is exactly how
these definitions drift apart.

Dependencies: base R and `stats`. Nothing else, deliberately.

## Provenance

The five working functions are a port of `R/spectral_dematel.R` from the
research code behind the source paper, which is the single source of truth for
every metric definition. The port is verbatim; everything added here is
documentation, tests and packaging.

The test suite compares each definition against an **independent route** to the
same number: values from the paper's precomputed corpus workbook, closed-form
spectra for matrices whose eigenvalues can be worked out on paper, the package's
own two input branches checked against each other, and finite differences for
the closed-form derivative. Agreement between a function and a literal it
produced proves only that nothing changed.

See [DEVELOPMENT.md](DEVELOPMENT.md) before changing anything in `R/`.

## License

MIT © Pavlos Delias
