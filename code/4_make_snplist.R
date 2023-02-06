#### Gather LD pruned SNPs
library(tidyverse)
library(data.table)
options("scipen"=100, "digits"=4)

args <- commandArgs(trailingOnly=TRUE)
datname <- args[1]
pruned_dat <- args[2]

setwd(paste0("/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/3_ldprune/", pruned_dat, "/pruned"))
allvar <- fread(paste0("/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/2_varqc/", datname, ".gz"))

#### After LD pruning, gather retained list ###

prunedall <- tibble()
for (idx in 1:22){
  pruned <- fread(paste0("chr", idx, ".pruned"))
  prunedall <- prunedall %>% bind_rows(pruned)
  print(paste("gather pruned out list from chr:", idx))
}
prunedall <- prunedall %>% 
  separate(col="variant", into=c("chr", "pos", "ref", "alt"), sep=":") %>% 
  mutate(chr = as.numeric(chr), pos = as.numeric(pos))
nrow(prunedall)

prunedall <- prunedall %>% inner_join(allvar, by=c("chr", "pos", "ref", "alt"))
prunedall %>% distinct(chr, pos, ref, alt)
prunedall %>% distinct(rsid)

prunedall %>% write_tsv("./allvar.pruned")
