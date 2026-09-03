#!/usr/bin/env Rscript
#
# Times TreeDist's MutualClusteringInfo and ClusteringInfoDistance on the trees written by
# run.jl. Both sides read the same Newick files and exclude parsing, so the benchmark
# compares only the tree computations.

suppressMessages(library(TreeDist))
suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
treedir <- args[[1]]
outfile <- args[[2]]

# Repeat a call until the sample spans the clock's resolution.
time_median <- function(f, budget = 0.5, min_reps = 5L) {
  times <- numeric(0)
  start <- Sys.time()
  repeat {
    t0 <- Sys.time()
    f()
    times <- c(times, as.numeric(Sys.time() - t0, units = "secs"))
    if (length(times) >= min_reps &&
        as.numeric(Sys.time() - start, units = "secs") > budget) break
    if (length(times) >= 1000L) break
  }
  median(times)
}

# R exposes peak memory through its garbage collector. This is not comparable to Julia's
# allocation count.
gc_megabytes <- function(f) {
  invisible(gc(reset = TRUE, full = TRUE))
  f()
  sum(gc(full = TRUE)[, "max used"] * c(8, 8)) / 1024^2
}

rows <- list()

for (n in c(10L, 50L, 200L, 1000L)) {
  f1 <- file.path(treedir, sprintf("pair_%d_a.nwk", n))
  if (!file.exists(f1)) next
  t1 <- read.tree(f1)
  t2 <- read.tree(file.path(treedir, sprintf("pair_%d_b.nwk", n)))

  fns <- list(
    mci = function() MutualClusteringInfo(t1, t2),
    cid = function() ClusteringInfoDistance(t1, t2)
  )
  for (metric in names(fns)) {
    f <- fns[[metric]]
    message(metric, " at ", n, " taxa")
    rows[[length(rows) + 1L]] <- data.frame(
      metric = metric, case = "pair", n = n, ntrees = 2L,
      seconds = time_median(f), megabytes = gc_megabytes(f), value = as.numeric(f())
    )
  }
}

write.table(do.call(rbind, rows), outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("TreeDist", as.character(packageVersion("TreeDist")),
    "| ape", as.character(packageVersion("ape")),
    "| R", paste0(R.version$major, ".", R.version$minor), "\n")
