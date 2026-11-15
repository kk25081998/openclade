# Cross-check OpenClade against TreeSearch (Smith 2023; Morphy library by
# Brazeau) on the 31 published matrices distributed with TreeSearch.
#
# For each matrix, OpenClade's best tree and 25 random trees are scored with
# TreeSearch::TreeLength() and compared with the scores OpenClade reported, for
# three analyses, all with the inapplicable algorithm of Brazeau, Guillerme &
# Smith (2019):
#   1. equal weights;
#   2. implied weights, concavity k = 10;
#   3. equal weights with some multistate characters ordered (TreeTools::Decompose()),
#      none of them scored with the inapplicable algorithm;
#   4. as 3, but including characters scored with the inapplicable algorithm. Here
#      OpenClade's default counts each extra applicable region once per character,
#      whereas a decomposed character counts it once per binary column, so this
#      comparison uses OpenClade's TreeSearch-compatible mode (see VALIDATION.md).
#
# Run in R or RStudio:
#   source("https://raw.githubusercontent.com/kk25081998/openclade/main/validation/treesearch/check.R")

base <- "https://raw.githubusercontent.com/kk25081998/openclade/main/validation/treesearch/"
treesearchCommit <- "78ee3d969fb7ab79eb054789dc4ec37da331ba3e" # source of the .nex files used by OpenClade

for (pkg in c("TreeTools", "TreeSearch")) {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    message("Installing ", pkg, " from CRAN...")
    install.packages(pkg, repos = "https://cloud.r-project.org")
  }
}
suppressPackageStartupMessages({
  library("TreeTools")
  library("TreeSearch")
})

# Use local files when run from a checkout (e.g. GitHub Actions), else fetch them.
useLocal <- file.exists("expected.csv") && dir.exists("trees")
src <- function(f) if (useLocal) f else url(paste0(base, f))
expected <- read.csv(src("expected.csv"), stringsAsFactors = FALSE)
datasets <- unique(expected$dataset)
results <- list()
norm <- function(x) tolower(gsub("[^A-Za-z0-9]", "", x))

for (ds in datasets) {
  nexFile <- system.file("datasets", paste0(ds, ".nex"), package = "TreeSearch")
  if (!nzchar(nexFile)) {
    nexFile <- tempfile(fileext = ".nex")
    download.file(paste0("https://raw.githubusercontent.com/ms609/TreeSearch/",
                         treesearchCommit, "/inst/datasets/", ds, ".nex"),
                  nexFile, quiet = TRUE)
  }
  dat <- ReadAsPhyDat(nexFile)
  ocTaxa <- readLines(src(paste0("trees/", ds, ".taxa.txt")))
  # Map OpenClade's row labels (t1..tN, file order) to TreeSearch's taxon names
  idx <- match(norm(ocTaxa), norm(names(dat)))
  if (anyNA(idx)) stop(ds, ": could not match taxa ", paste(ocTaxa[is.na(idx)], collapse = ", "))
  trees <- ape::read.tree(text = readLines(src(paste0("trees/", ds, ".tre"))))
  trees <- lapply(trees, function(tr) {
    tr$tip.label <- names(dat)[idx[as.integer(sub("^t", "", tr$tip.label))]]
    tr
  })
  class(trees) <- "multiPhylo"
  ts <- TreeLength(trees, dat)
  ex <- expected[expected$dataset == ds, ]
  # implied weights: fits agree to the 6 decimals OpenClade reports
  iw <- TreeLength(trees, dat, concavity = 10)
  # ordered characters, numbered as in the matrix: line 1 = set A, line 2 = set B
  ordLines <- readLines(src(paste0("trees/", ds, ".ordered.txt")))
  sets <- lapply(ordLines[1:2], function(l) as.integer(strsplit(trimws(l), " +")[[1]]))
  sets <- lapply(sets, function(x) x[!is.na(x)])
  ordA <- if (length(sets[[1]])) TreeLength(trees, Decompose(dat, sets[[1]])) else ts
  ordB <- if (length(sets[[2]])) TreeLength(trees, Decompose(dat, sets[[2]])) else ts
  results[[ds]] <- data.frame(dataset = ds, tree = ex$tree, kind = ex$kind,
                              openclade = ex$openclade_length, treesearch = ts,
                              match = ex$openclade_length == ts,
                              openclade_iw10 = ex$openclade_iw10, treesearch_iw10 = round(iw, 6),
                              match_iw10 = abs(ex$openclade_iw10 - iw) < 2e-6,
                              openclade_ordered = ex$openclade_ordered, treesearch_ordered = ordA,
                              match_ordered = ex$openclade_ordered == ordA,
                              openclade_ordered_na = ex$openclade_ordered_na, treesearch_ordered_na = ordB,
                              match_ordered_na = ex$openclade_ordered_na == ordB)
  cat(sprintf("%-14s equal weights %2d/%2d   implied weights %2d/%2d   ordered (%2d chars) %2d/%2d   ordered incl. inapplicable (%2d chars) %2d/%2d\n",
              ds, sum(results[[ds]]$match), nrow(ex), sum(results[[ds]]$match_iw10), nrow(ex),
              length(sets[[1]]), sum(results[[ds]]$match_ordered), nrow(ex),
              length(sets[[2]]), sum(results[[ds]]$match_ordered_na), nrow(ex)))
}

all <- do.call(rbind, results)
outFile <- file.path(getwd(), "openclade_vs_treesearch.csv")
write.csv(all, outFile, row.names = FALSE)
cat("\n", sum(all$match), " of ", nrow(all), " trees: OpenClade length == TreeSearch length (equal weights)\n", sep = "")
cat(sum(all$match_iw10), " of ", nrow(all), " trees: OpenClade fit == TreeSearch fit (implied weights, k = 10)\n", sep = "")
cat(sum(all$match_ordered), " of ", nrow(all), " trees: OpenClade length == TreeSearch length (ordered characters)\n", sep = "")
cat(sum(all$match_ordered_na), " of ", nrow(all), " trees: OpenClade (TreeSearch-compatible mode) == TreeSearch (ordered characters incl. inapplicable)\n", sep = "")
cat("TreeSearch ", format(packageVersion("TreeSearch")), ", TreeTools ", format(packageVersion("TreeTools")),
    ", ", R.version.string, "\n", sep = "")
cat("Full table written to openclade_vs_treesearch.csv\n")
nBad <- sum(!all$match) + sum(!all$match_iw10) + sum(!all$match_ordered) + sum(!all$match_ordered_na)
if (nBad > 0) {
  cat("MISMATCHES: ", nBad, "\n", sep = "")
  quit(status = 1)
}
