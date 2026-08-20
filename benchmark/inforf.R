#!/usr/bin/env Rscript
#
# Times TreeDist's InfoRobinsonFoulds on the trees written by run.jl, and appends one row
# per case to results_inforf_r.tsv. Both sides read the same Newick files and both exclude
# parsing, so only the distance computation is compared.

suppressMessages(library(TreeDist))
suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
treedir <- args[[1]]
outfile <- args[[2]]

# Repeats the call until at least `budget` seconds have elapsed, so that operations far
# faster than the clock's resolution are still timed meaningfully.
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

# R reports memory only through the garbage collector, so this is the peak the GC saw
# rather than a total allocation count. It is not comparable to Julia's figure.
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

  secs <- time_median(function() InfoRobinsonFoulds(t1, t2))
  mb <- gc_megabytes(function() InfoRobinsonFoulds(t1, t2))
  value <- InfoRobinsonFoulds(t1, t2)
  rows[[length(rows) + 1L]] <- data.frame(
    case = "pair", n = n, ntrees = 2L, seconds = secs, megabytes = mb, value = value
  )
}

write.table(do.call(rbind, rows), outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("TreeDist", as.character(packageVersion("TreeDist")),
    "| ape", as.character(packageVersion("ape")),
    "| R", paste0(R.version$major, ".", R.version$minor), "\n")
