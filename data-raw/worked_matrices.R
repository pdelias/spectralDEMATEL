# Builds data/worked_matrices.rda from the transcribed CSVs in ../reference/.
#
# Run from the package root:  Rscript data-raw/worked_matrices.R
#
# The CSVs are headerless, comma separated, one row per factor. They come from
# the JDS paper repository's data/raw/matrices/ and are copied into
# ../reference/matrices/ so this project does not depend on that repository.

REF <- "../reference/matrices"

read_matrix <- function(file, doi, label) {
  A <- as.matrix(utils::read.csv(file.path(REF, file), header = FALSE))
  storage.mode(A) <- "double"
  dimnames(A) <- NULL
  attr(A, "doi") <- doi
  attr(A, "label") <- label
  A
}

worked_matrices <- list(
  fefo_stock_control = read_matrix(
    "10.1016_j.jafr.2025.101848.csv",
    doi   = "10.1016/j.jafr.2025.101848",
    label = "FEFO stock control"),
  resilience_capabilities = read_matrix(
    "10.1007_s12063-024-00470-8.csv",
    doi   = "10.1007/s12063-024-00470-8",
    label = "Resilience capabilities")
)

stopifnot(
  dim(worked_matrices$fefo_stock_control) == c(13, 13),
  dim(worked_matrices$resilience_capabilities) == c(15, 15)
)

dir.create("data", showWarnings = FALSE)
save(worked_matrices, file = "data/worked_matrices.rda",
     version = 2, compress = "xz")

cat("wrote data/worked_matrices.rda\n")
