#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(dada2)

fwd_in <- snakemake@input[["r1"]]
rev_in <- snakemake@input[["r2"]]

fwd_out <- snakemake@output[["r1"]]
rev_out <- snakemake@output[["r2"]]
 
# debugging
# fwd_in <- "data/reads/J2HBL-3633-248-00-01_S72_L001_R1_001.fastq.gz"
# rev_in <- "data/reads/J2HBL-3633-248-00-01_S72_L001_R2_001.fastq.gz"
# fwd_out <- "test1.fastq.gz"
# rev_out <- "test2.fastq.gz"
# data.table::fread(paste0("zcat ", fwd_in, " | head"))
# data.table::fread(paste0("zcat ", rev_in, " | head"))

filterAndTrim(fwd = fwd_in,
              filt = fwd_out,
              rev = rev_in,
              filt.rev = rev_out,
              maxN = 0,
              multithread = FALSE,
              matchIDs = TRUE)

sessionInfo()
