#!/usr/bin/env Rscript
#
# Times Quartet's quartet distance on the trees written by run.jl, and appends one row per
# case to the output TSV. Both sides read the same Newick files and both exclude parsing, so
# only the distance computation is compared.
#
# The distance itself is reported alongside the timing so that run.jl can check the two
# implementations agree on every tree pair it benchmarked.

suppressMessages(library(Quartet))
suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
treedir <- args[[1]]
outfile <- args[[2]]
sizes <- as.integer(strsplit(args[[3]], ",")[[1]])

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

# A quartet resolved by only one tree counts as a difference, so the distance is the
# conflicting quartets plus those resolved on one side alone.
#
# The `N` column is not read: it is 2 * Q, which overflows R's 32-bit integers at 477 tips
# while Q itself still fits. The warning that overflow raises is expected and suppressed.
quartet_distance <- function(t1, t2) {
  st <- suppressWarnings(QuartetStatus(c(t1, t2))[2, ])
  as.numeric(st[["d"]]) + as.numeric(st[["r1"]]) + as.numeric(st[["r2"]])
}

rows <- list()

for (n in sizes) {
  f1 <- file.path(treedir, sprintf("pair_%d_a.nwk", n))
  if (!file.exists(f1)) next
  t1 <- read.tree(f1)
  t2 <- read.tree(file.path(treedir, sprintf("pair_%d_b.nwk", n)))

  value <- quartet_distance(t1, t2)
  secs <- time_median(function() quartet_distance(t1, t2))
  mb <- gc_megabytes(function() quartet_distance(t1, t2))
  rows[[length(rows) + 1L]] <- data.frame(
    case = "quartet", n = n, ntrees = 2L, seconds = secs, megabytes = mb, value = value
  )
}

write.table(do.call(rbind, rows), outfile, sep = "\t", row.names = FALSE, quote = FALSE)
cat("Quartet", as.character(packageVersion("Quartet")),
    "| ape", as.character(packageVersion("ape")),
    "| R", paste0(R.version$major, ".", R.version$minor), "\n")
