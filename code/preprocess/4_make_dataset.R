### Gather summary statistics from files
library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
filename <- args[1] # phe file
datname <- args[2] # sig_pval5e-08_maf1e-02_nct2_1kg
querylist <- args[3] # allvar.pruned

indir <- "./PanUKBB/"
outpath <- paste0("./data/PanUKBB/4_sumstats/", datname, "/")

## query file need chr, pos, ref, alt, rsid
querylist_path <- paste0("./data/PanUKBB/3_ldprune/", datname, "/pruned/", querylist)
querydf <- fread(querylist_path) %>% mutate(chr=as.numeric(chr), pos=as.numeric(pos))

## test file: continuous-22409-both_sexes-irnt.tsv.bgz
dat <- read_delim(gzfile(paste0(indir, filename)), delim='\t') %>%
  mutate(chr=as.numeric(chr), pos=as.numeric(pos))

# some contains Inf Z score b/c SE = 0 (will force these to be NA)
dat_joined <- inner_join(querydf, dat, by=c("chr", "pos", "ref", "alt")) %>%
  drop_na(beta_EUR, se_EUR) %>%
  mutate(z_EUR = beta_EUR/se_EUR) %>%
  select(chr, pos, ref, alt, rsid, varid, z_EUR)

dat_joined %>% fwrite(paste0(outpath, filename), sep='\t')
