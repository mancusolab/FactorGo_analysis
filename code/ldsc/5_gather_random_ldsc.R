library(tidyverse)
library(data.table)
options("scipen"=100, "digits"=4)

resdir <- paste0("./results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/enrich/random_gene_annotations/ldsc_cts/")
outdir <- "./results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/enrich/random_gene_annotations/"

num_fac <- 100

######## gather all result in one file #######

vb_allres <- tibble()
svd_allres <- tibble()

# HERE use wt_cosine of VB results
for (fac in 1:num_fac){
  vbfile <- fread(paste0(resdir,
                         "vb.scaleTrue.K100.F", fac, ".pseudoNTRUE.wt_cosine.cell_type_results.txt"),
                  header=TRUE)
  vb_allres <- bind_rows(vb_allres, vbfile)

  svdfile <- fread(paste0(resdir,
                          "svd.scalesnps.K100.F", fac, ".pseudoNTRUE.cell_type_results.txt"),
                   header=TRUE)
  svd_allres <- bind_rows(svd_allres, svdfile)
}

vb_allres %>%
  mutate(fac = rep(1:num_fac, each=10)) %>%
  write_tsv(paste0(outdir,"vb.scaleTrue.K", num_fac, ".randomgene.pseudoNTRUE.wt_cosine.4400blocks.cell_type_results.txt"))

svd_allres %>%
  mutate(fac = rep(1:num_fac, each=10)) %>%
  write_tsv(paste0(outdir, "svd.scalesnps.K", num_fac, ".randomgene.pseudoNTRUE.4400blocks.cell_type_results.txt"))
