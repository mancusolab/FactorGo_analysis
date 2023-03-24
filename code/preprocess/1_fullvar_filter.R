#### Variant filter for full variant file ####
## data.table tutorial: https://atrebas.github.io/post/2019-03-03-datatable-dplyr/
library(data.table)
library(tidyverse)
options("scipen"=100, "digits"=4)

args <- commandArgs(trailingOnly=TRUE)
info_pass <- as.numeric(args[1])
maf <- as.numeric(args[2])
high_qual <- as.logical(args[3])
outpath <- args[4]

fullvarDT <- fread("./PanUKBB/full_variant_qc_metrics.tsv")
names(fullvarDT)[1] <- "chr"

# filter by autosomal variant
fullvarDT <- fullvarDT[chr %chin% as.character(1:22)] # fast %in% for character
print(table(fullvarDT[, chr]))
print(paste0("Subset to ", nrow(fullvarDT), " autosomal variants"))

# filter by info score, maf, high quality variants
fullvarDT[, maf_EUR := ifelse(af_EUR > 0.5, 1-af_EUR, af_EUR)]
summary(fullvarDT[, maf_EUR])

DT_filtered <- fullvarDT[(info > info_pass & maf_EUR > maf & high_quality == high_qual)]
summary(DT_filtered[, maf_EUR])
print(paste0("Subset to ", nrow(DT_filtered), "  variants with maf >", maf, " and high_qual:", high_qual))

DT_filtered %>% fwrite(paste0(outpath, "full_variant_info", info_pass, "_maf", maf, "_highqual", high_qual, ".gz"), sep='\t')
