#!/usr/bin/env Rscript

log <- file(snakemake@log[[1]], open = "wt")
sink(log, type = "message")
sink(log, type = "output", append = TRUE)

library(dada2)

r1_files <- snakemake@input[['r1']]
r2_files <- snakemake@input[['r2']]
ref_fasta <- snakemake@input[['fa']]

threads <- snakemake@threads[[1]]

seqtab_file <- snakemake@output[["seqtab"]]
nochim_file <- snakemake@output[["nochim"]]
taxa_file <- snakemake@output[["taxa"]]
errF_plot_file <- snakemake@output[["errF_plot"]]
errR_plot_file <- snakemake@output[["errR_plot"]]

# dev
# r1_files <- list.files("output/020_preprocess",
#            pattern = "cutadapt_r1",
#            full.names = TRUE)[1:4]
# r2_files <- list.files("output/020_preprocess",
#                        pattern = "cutadapt_r2",
#                        full.names = TRUE)[1:4]
# threads <- 8
# ref_fasta <- "data/its1-58024_dada2.fa"

# set up trim files
names(r1_files) <- sub("_cutadapt_r1.fq.gz", "", basename(r1_files))
names(r2_files) <- sub("_cutadapt_r2.fq.gz", "", basename(r2_files))

mytmp <- function(x, fileext) {tempfile(fileext = fileext)}
r1_tmp <- sapply(r1_files, mytmp, fileext = ".fq.gz")
r2_tmp <- sapply(r2_files, mytmp, fileext = ".fq.gz")

# rerun filter... not sure why?
out <- filterAndTrim(
    r1_files,
    r1_tmp,
    r2_files,
    r2_tmp,
    maxN = 0,
    maxEE = c(2, 2), 
    truncQ = 2,
    minLen = 50,
    rm.phix = TRUE,
    compress = TRUE,
    multithread = threads,
    matchIDs = TRUE) 

# learn errors
errF <- learnErrors(r1_tmp, multithread = threads)
errR <- learnErrors(r2_tmp, multithread = threads)
errF_plot <- plotErrors(errF, nominalQ = TRUE)
errR_plot <- plotErrors(errR, nominalQ = TRUE)

# dereplicate
derepF <- derepFastq(r1_tmp[file.exists(r1_tmp)], verbose = TRUE)
derepR <- derepFastq(r2_tmp[file.exists(r2_tmp)], verbose = TRUE)

# run dada algo
dadaF <- dada(derepF, err = errF, multithread = threads)
dadaR <- dada(derepR, err = errR, multithread = threads)

# merge reads
mergers <- mergePairs(dadaF,
                      derepF,
                      dadaR,
                      derepR,
                      verbose = TRUE,
                      justConcatenate = TRUE)

# make output
seqtab <- makeSequenceTable(mergers)

# remove chimeras
seqtab.nochim <- removeBimeraDenovo(seqtab,
                                    method = "consensus",
                                    multithread = threads,
                                    verbose = TRUE)

# assign taxonomy
taxa <- assignTaxonomy(seqtab.nochim,
                       ref_fasta,
                       multithread = threads)

# write output
saveRDS(seqtab, seqtab_file)
saveRDS(seqtab.nochim, nochim_file)
saveRDS(taxa, taxa_file)
ggplot2::ggsave(errF_plot_file,
                errF_plot,
                width = 10,
                height = 7.5,
                units = "in",
                device = cairo_pdf)
ggplot2::ggsave(errR_plot_file,
                errR_plot,
                width = 10,
                height = 7.5,
                units = "in",
                device = cairo_pdf)

# log
sessionInfo()
