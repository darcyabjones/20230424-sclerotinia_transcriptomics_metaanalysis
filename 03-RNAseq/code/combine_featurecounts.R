#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)

indir <- args[1] # "output/feature_counts/"
outfile <- args[2]

if (length(args) < 2) {
  cat("USAGE: script.R indir outfile\n", file=stderr())
  quit(save = "no", status = 1)
} 

l <- list.files(path = indir, pattern = "*.txt$")
files <- list()

for (f in l) {
  t <- read.table(paste0(indir, "/", f), header = 2);
  n <- gsub(f, pattern = ".txt", replacement = "");
  colnames(t)[ncol(t)] <- n;
  files[[n]] <- t
}
out <- Reduce(f = function(x, y) {merge(x, y, by = c("Geneid", "Chr", "Start", "End", "Strand", "Length"))}, x = files)

write.table(out, outfile, quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)


l <- list.files(path = indir, pattern = "*.txt.summary$")
files <- list()

for (f in l) {
  t <- read.table(paste0(indir, "/", f), header = 1);
  n <- gsub(f, pattern = ".txt.summary", replacement = "");
  colnames(t)[ncol(t)] <- n;
  files[[n]] <- t
}
out <- Reduce(f = function(x, y) {merge(x, y, by = "Status")}, x = files)

write.table(out, paste0(outfile, ".summary"), quote = FALSE, sep="\t", row.names = FALSE, col.names = TRUE)
