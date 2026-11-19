# Builds R/corpus-h.R, the aggregate constants of the 81-model corpus that the
# H reading is read against.
#
# Run from the package root:
#
#   Rscript data-raw/corpus_constants.R [objects] [library]
#
#   objects  folder holding results.rds, quadrants.rds and corpus.rds, copies
#            of the targets `results`, `quadrants` and `corpus` of the paper's
#            pipeline. Default: ../reference/corpus/revised-2026-11-01
#   library  an R library that holds quantreg, which the pipeline uses for the
#            trade-off line. The paper repository's renv library has it.
#            Without this argument quantreg must be installed.
#
# No constant below is typed by hand. Each one is either a key of the paper's
# results or is computed from the pipeline's per-model table and checked
# against a key. The script stops if any check fails.
#
# Only aggregates leave this script. The per-model rows stay in ../reference/,
# outside the package.

args <- commandArgs(TRUE)
objects <- if (length(args) >= 1) args[1] else "../reference/corpus/revised-2026-11-01"
if (length(args) >= 2) .libPaths(c(args[2], .libPaths()))
if (!requireNamespace("quantreg", quietly = TRUE)) {
  stop("quantreg is needed to refit the trade-off line; give a library that holds it")
}

res  <- readRDS(file.path(objects, "results.rds"))
quad <- readRDS(file.path(objects, "quadrants.rds"))
corp <- readRDS(file.path(objects, "corpus.rds"))
d <- corp$diagnostics

types <- c("diffuse-amplified", "hierarchical-dampened",
           "hierarchical-amplified", "diffuse-dampened")
short <- c(`diffuse-amplified` = "da", `hierarchical-dampened` = "hd",
           `hierarchical-amplified` = "ha", `diffuse-dampened` = "dd")

# ---- the run these objects come from -----------------------------------------

stopifnot(
  identical(res$sample, "eligible"),
  identical(corp$sample, "eligible"),
  identical(res$quad_axis, "H"),
  !isTRUE(res$provisional),
  res$htype_cut == 0.10,
  res$quad_mu_cut == 0.5,
  nrow(d) == res$htype_n,
  all(is.finite(d$mu_max)), all(is.finite(d$hierarchy_H)),
  all(is.finite(d$dominance))
)

# ---- counts and multipliers: keys of the results -----------------------------

count <- vapply(types, function(ty) {
  as.integer(res[[paste0("htype_", short[[ty]], "_n")]])
}, integer(1))
multiplier <- vapply(types, function(ty) {
  res[[paste0("quad_", short[[ty]], "_multiplier")]]
}, numeric(1))

cells <- quad$cells
stopifnot(
  sum(count) == res$htype_n,
  all(count == cells$n[match(types, cells$quadrant)]),
  all(abs(multiplier -
            round(cells$median_multiplier[match(types, cells$quadrant)], 2)) < 1e-12)
)

# The types recomputed from the per-model table, with the pipeline's rule
# (strictly above each cut).
type_of <- function(mu, h) {
  paste0(ifelse(h > 0.10, "hierarchical", "diffuse"), "-",
         ifelse(mu > 0.5, "amplified", "dampened"))
}
observed_type <- type_of(d$mu_max, d$hierarchy_H)
stopifnot(all(table(observed_type)[types] == count))

# ---- quantiles: computed from the per-model table ------------------------------

# The pipeline stores the quartiles of H (quadrants$hier_q) and three medians.
# The other percentiles are not keys, so they are computed here with the
# quantile type the pipeline uses throughout (R's default, type 7).
probs <- c(0.05, 0.25, 0.50, 0.75, 0.95)
qs <- function(x) unname(stats::quantile(x, probs))
quantiles <- list(mu_max      = qs(d$mu_max),
                  hierarchy_h = qs(d$hierarchy_H),
                  dominance   = qs(d$dominance))
stopifnot(
  max(abs(quantiles$hierarchy_h[2:4] - quad$hier_q)) < 1e-12,
  round(quantiles$mu_max[3], 3)      == res$median_mu_max,
  round(quantiles$hierarchy_h[3], 3) == res$hier_H_median,
  round(quantiles$dominance[3], 3)   == res$dominance_median,
  round(quantiles$hierarchy_h[4], 3) == res$quad_hier_q75
)
quantiles <- lapply(quantiles, round, 4)

# ---- the trade-off line: refitted, then checked against four keys ---------------

# tradeoff_reconciliation() in the pipeline's task11_quadrants.R fits a median
# regression of H on coupling and keeps the fitted values only. The same call
# is repeated here to recover the intercept and the slope.
fit <- quantreg::rq(hierarchy_H ~ mu_max, tau = 0.5, data = d)
hat <- as.numeric(stats::predict(fit, d))
resid_sd <- stats::sd(d$hierarchy_H - hat)
scal <- quad$reconcile$scalars

expected_type <- type_of(d$mu_max, hat)
wrong <- expected_type != observed_type

stopifnot(
  identical(scal$method, "median regression"),
  round(unname(stats::coef(fit)[["mu_max"]]), 3) == res$tradeoff_q50_slope,
  abs(resid_sd - scal$resid_sd) < 1e-10,
  sum(wrong) == res$quad_alone_wrong_n,
  sum(wrong & observed_type == "hierarchical-amplified") == res$quad_alone_wrong_ha_n,
  sum(wrong & observed_type == "diffuse-dampened") == res$quad_alone_wrong_dd_n,
  round(scal$r2, 2) == res$quad_tradeoff_r2
)

tradeoff <- list(
  intercept = round(unname(stats::coef(fit)[["(Intercept)"]]), 6),
  slope     = round(unname(stats::coef(fit)[["mu_max"]]), 6),
  resid_sd  = round(resid_sd, 6),
  # Four decimals, so that rounding again to two gives the 0.75 the manuscript
  # prints. Three would store 0.745, which rounds down.
  r2        = round(scal$r2, 4),
  spearman  = res$tradeoff_rho_H,
  method    = scal$method
)

stopifnot(round(tradeoff$r2, 2) == res$quad_tradeoff_r2)

# The rounded line must still give the counts the manuscript reports.
hat_r <- tradeoff$intercept + tradeoff$slope * d$mu_max
stopifnot(sum(type_of(d$mu_max, hat_r) != observed_type) == res$quad_alone_wrong_n)

# ---- assemble and write ----------------------------------------------------------

CORPUS_H <- list(
  n          = as.integer(res$htype_n),
  tradeoff   = tradeoff,
  count      = count,
  multiplier = multiplier,
  quantiles  = quantiles,
  probs      = probs,
  # What the trade-off gets wrong when it is used alone to assign types.
  alone = list(wrong = as.integer(res$quad_alone_wrong_n),
               hierarchical_amplified = as.integer(res$quad_alone_wrong_ha_n),
               diffuse_dampened = as.integer(res$quad_alone_wrong_dd_n)),
  # Share of the attainable value that uniform effort captures, 1 - H: its
  # median, and the percentage of models that keep the loss within the cut.
  uniform_share = list(median = res$uniform_share_median,
                       within_cut_pct = res$uniform_share_ge90_pct),
  source = list(generated_at = res$generated_at,
                sample = res$sample,
                results_md5 = unname(tools::md5sum(file.path(objects, "results.rds"))))
)

header <- c(
  "# Generated by data-raw/corpus_constants.R. Do not edit by hand.",
  "#",
  "# Aggregate constants of the 81-model corpus of the revised source paper, on",
  "# which the H reading is cut. Each value is a key of the paper's results or is",
  "# computed from its per-model table and checked against a key. The script",
  "# documents the source of every one.",
  "",
  "#' Aggregate constants of the 81-model corpus, for the H reading.",
  "#' @noRd"
)
# Every value was rounded above, so fifteen significant digits print it exactly.
body <- deparse(CORPUS_H, width.cutoff = 60L,
                control = c("keepInteger", "niceNames"))
body[1] <- paste("CORPUS_H <-", body[1])
writeLines(c(header, body), "R/corpus-h.R")
stopifnot(identical(local({ source("R/corpus-h.R", local = TRUE); CORPUS_H }),
                    CORPUS_H))
cat("wrote R/corpus-h.R\n")
str(CORPUS_H)
