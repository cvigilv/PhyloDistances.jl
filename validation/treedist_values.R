#!/usr/bin/env Rscript
#
# Computes TreeDist's values for the tree pairs written by crosscheck.jl.
#
# Values are written with %.17g and compared on the Julia side. R's own as.numeric does not
# reliably round-trip that output, so parsing must not happen here.

suppressMessages(library(TreeDist))
suppressMessages(library(ape))

args <- commandArgs(trailingOnly = TRUE)
cases <- read.delim(args[[1]], stringsAsFactors = FALSE, colClasses = "character")
out <- file(args[[2]], "w")

cat("label\tnw1\tnw2\trf\tnorm\n", file = out)
for (i in seq_len(nrow(cases))) {
  t1 <- read.tree(text = cases$nw1[i])
  t2 <- read.tree(text = cases$nw2[i])
  rf <- RobinsonFoulds(t1, t2)
  nm <- RobinsonFoulds(t1, t2, normalize = TRUE)
  cat(cases$label[i], cases$nw1[i], cases$nw2[i],
      sprintf("%d", as.integer(rf)),
      if (is.nan(nm)) "NaN" else sprintf("%.17g", nm),
      sep = "\t", file = out)
  cat("\n", file = out)
}
close(out)

cat("TreeDist ", as.character(packageVersion("TreeDist")),
    ", ape ", as.character(packageVersion("ape")),
    ", R ", paste0(R.version$major, ".", R.version$minor), "\n", sep = "")
