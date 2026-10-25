# Cross-check OpenClade against TreeSearch (Smith 2023; Morphy library by
# Brazeau) on the 31 published matrices distributed with TreeSearch.
#
# For each matrix, OpenClade's best tree and 25 random trees are scored with
# TreeSearch::TreeLength() (equal weights, inapplicable algorithm of Brazeau,
# Guillerme & Smith 2019) and compared with the lengths OpenClade reported.
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
  results[[ds]] <- data.frame(dataset = ds, tree = ex$tree, kind = ex$kind,
                              openclade = ex$openclade_length, treesearch = ts,
                              match = ex$openclade_length == ts)
  cat(sprintf("%-14s %2d/%2d trees match   (best tree: OpenClade %d, TreeSearch %d)\n",
              ds, sum(results[[ds]]$match), nrow(ex),
              ex$openclade_length[ex$kind == "best"], ts[ex$kind == "best"]))
}

all <- do.call(rbind, results)
outFile <- file.path(getwd(), "openclade_vs_treesearch.csv")
write.csv(all, outFile, row.names = FALSE)
cat("\n", sum(all$match), " of ", nrow(all), " trees: OpenClade length == TreeSearch length\n", sep = "")
cat("TreeSearch ", format(packageVersion("TreeSearch")), ", TreeTools ", format(packageVersion("TreeTools")),
    ", ", R.version.string, "\n", sep = "")
cat("Full table written to ", outFile, "\n", sep = "")
