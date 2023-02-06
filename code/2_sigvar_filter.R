#### Combine list of sig variant from each phe file ###

library(data.table)
library(tidyverse)

args <- commandArgs(trailingOnly=TRUE)
pval <- as.numeric(args[1])
nct <- as.numeric(args[2])
maf <- as.numeric(args[3])
phe_file <- args[4]
options("scipen"=100, "digits"=4)

outpath <- "/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/2_varqc/"
sigdir <- "/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/1_sigSNP/"
filenames <- read_tsv(phe_file) %>% pull(filename)

## use QCed full variant data
fullvarqc <- fread("/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/2_varqc/full_variant_info0.9_maf1e-04_highqualTRUE.gz")

## if need further prune by maf
fullvarqc <- fullvarqc %>% 
  filter(maf_EUR > maf) %>% 
  mutate(chr = as.numeric(chr))

### combine significant variants (no filtering on INFO, MAF, high_qual)
alldat <- tibble()
for (idx in seq_along(filenames)){
  dat <- fread(paste0(sigdir, filenames[idx])) %>% 
    filter(chr %in% 1:22) %>% 
    mutate(chr = as.numeric(chr)) %>% 
    filter(pval_EUR < log(pval))
  
  if (nrow(dat) > 0){
    alldat <- bind_rows(alldat, dat)
  }
  if (idx %% 100 == 0){
    print(paste("Read in", idx))
  }
}

## restrict to #associations threshold and filter by INFO, MAF, high_qual
alldat %>% 
  group_by(chr, pos, ref, alt) %>% 
  add_count() %>% 
  filter(n >= nct) %>%
  ungroup() %>% 
  inner_join(fullvarqc, by=c("chr", "pos", "ref", "alt")) %>% 
  distinct(chr, pos, ref, alt, .keep_all = TRUE) %>% 
  select(chr, pos, ref, alt, rsid, varid, maf_EUR) %>% 
  fwrite(paste0(outpath, "sig", "_pval", pval, "_nct", nct, "_maf", maf, ".gz"), sep='\t')

