# Prints the reference values that the M1 tests hold as literals.
#
# Run from the package root, in a fresh R process:
#
#   Rscript data-raw/reference_values.R "<path to the paper repository>"
#
# Nothing here loads the development version of the package. The values come
# from two independent routes:
#
#   1. The paper's own code, R/baseline/spectral_baseline.R, sourced into its
#      own environment and run on the two transcribed matrices.
#   2. The released 0.6.0 of this package, as installed in the library. The
#      script stops if another version is installed.
#
# The output was pasted into tests/testthat/test-hierarchy-h.R, test-nulls.R
# and test-leverage.R on 2026-11-19. If a literal there has to change, rerun
# this script and say why in NEWS.md.

args <- commandArgs(TRUE)
if (length(args) < 1) stop("give the path to the paper repository")
jds <- args[1]
ref <- "../reference/matrices"

paper <- new.env()
sys.source(file.path(jds, "R/baseline/spectral_baseline.R"), paper)

read_matrix <- function(file) {
  A <- as.matrix(utils::read.csv(file.path(ref, file), header = FALSE))
  storage.mode(A) <- "double"
  dimnames(A) <- NULL
  A
}
fefo <- read_matrix("10.1016_j.jafr.2025.101848.csv")
resl <- read_matrix("10.1007_s12063-024-00470-8.csv")

show <- function(label, x) {
  cat("\n# ", label, "\n", sep = "")
  dput(x, control = "digits17")
}

# ---- route 1: the paper's code ---------------------------------------------

# The diagonal is kept, because the package keeps it. For these two matrices H
# is the same either way: the FEFO diagonal is constant and the resilience
# diagonal is zero.
for (nm in c("fefo", "resl")) {
  A <- get(nm)
  d <- paper$diagnostics(paper$dematel_normalise(A, diagonal = "keep"))
  show(paste(nm, "hierarchy_H, paper diagnostics()"), d[["hierarchy_H"]])
  show(paste(nm, "H_from_sd(), paper"),
       paper$H_from_sd(d[["hierarchy_sd"]], nrow(A)))
  z <- paper$diagnostics(paper$dematel_normalise(A, diagonal = "zero"))
  show(paste(nm, "hierarchy_H with the diagonal zeroed"), z[["hierarchy_H"]])
}

# The two nulls on the FEFO matrix, six draws each from a fixed seed. The
# package must reproduce these from the same seed, because it ports the same
# shuffles and draws them in the same order.
D <- paper$dematel_normalise(fefo, diagonal = "keep")
for (scheme in c("full", "rows")) {
  set.seed(20261119)
  nd <- paper$null_distribution(D, B = 6, null = scheme)
  show(paste("FEFO null,", scheme, ": coupling"), unname(nd[, "coupling"]))
  show(paste("FEFO null,", scheme, ": hierarchy_H"), unname(nd[, "hierarchy_H"]))
  show(paste("FEFO null,", scheme, ": spectral_ratio"),
       unname(nd[, "spectral_ratio"]))
}

# Mean, standard deviation and z of the paper's compare_to_null(), on the same
# six draws of the rows null.
set.seed(20261119)
nd <- paper$null_distribution(D, B = 6, null = "rows")
cmp <- paper$compare_to_null(paper$diagnostics(D), nd)
keep <- cmp$statistic %in% c("coupling", "hierarchy_H", "spectral_ratio",
                             "hierarchy_sd")
show("compare_to_null on the rows null: statistic, null_mean, null_sd, z",
     cmp[keep, c("statistic", "null_mean", "null_sd", "z")])

# Row and column totals of the total-relation matrix, for prominence and
# relation.
Tm <- paper$total_relation(D)
show("FEFO row totals of T (dispatch)", rowSums(Tm))
show("FEFO column totals of T (receipt)", colSums(Tm))

# ---- route 2: the released 0.6.0 --------------------------------------------

stopifnot(as.character(utils::packageVersion("spectralDEMATEL")) == "0.6.0")
ens <- spectralDEMATEL::surrogate_ensemble(fefo, B = 8, seed = 99)
show("0.6.0 surrogate_ensemble(fefo, B = 8, seed = 99)", as.list(ens))
show("its attributes: requested, attempts, complete",
     list(attr(ens, "requested"), attr(ens, "attempts"), attr(ens, "complete")))

sp <- spectralDEMATEL::surrogate_position(fefo, B = 25, seed = 99)
show("0.6.0 surrogate_position(fefo, B = 25, seed = 99)$metrics",
     as.list(sp$metrics))
show("its type, type_share, B, attempts",
     list(sp$type, sp$type_share, sp$B, sp$attempts))

ms <- spectralDEMATEL::measurement_stability(resl, B = 40, seed = 5)
show("0.6.0 measurement_stability(resl, B = 40, seed = 5): type and share",
     list(ms$observed_type, ms$share_same, ms$B_admissible))
