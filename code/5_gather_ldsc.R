library(tidyverse)
library(data.table)
options("scipen"=100, "digits"=4)

args <- commandArgs(trailingOnly=TRUE)
annot_name <- args[1] # SEG, "IMPACT707";
filelist <- args[2] # param_SEG_file_annot,  filelist="param_IMPACT_subset"
facnum <- args[3] # onefac, all
focal <- as.logical(args[4]) # True, False

resdir <- paste0("/project/nmancuso_8/elezhang/projects/FA/results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/enrich/", 
                 annot_name, "/ldsc_cts/")
outdir <- paste0("/project/nmancuso_8/elezhang/projects/FA/results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/enrich/", 
                 annot_name, "/")
filenames <- read_tsv(filelist, F) %>% distinct(X1) %>% pull(X1)

num_annot <- ifelse(annot_name == "SEG", 205, 707)
print(paste0("number of annotation: ", num_annot))

if (facnum == "onefac"){
  for (onefile in filenames){
    allres <- tibble()
    for (annot in 1:num_annot){
      file1 <- fread(paste0(resdir, onefile, ".Annot", annot, ".cell_type_results.txt"), header=T)
      allres <- bind_rows(allres, file1)
    } 
    allres %>%
      arrange(Coefficient_P_value) %>%
      write_tsv(paste0(resdir, onefile, ".4400blocks.cell_type_results.txt"))
    print(paste("Done with:", onefile))
  }
}
######## gather all result in one file #######

vb_allres <- tibble()
svd_allres <- tibble()

# HERE use wt_cosine of VB results
if (facnum == "all"){
  K <- 100
  for (fac in 1:K){
    vbfile <- fread(paste0(resdir,
                           "vb.scaleTrue.K100.F", fac, ".pseudoNTRUE.wt_cosine.4400blocks.cell_type_results.txt"),
                    header=TRUE) %>%
      mutate(Factor = fac)
    vb_allres <- bind_rows(vb_allres, vbfile)
    
    svdfile <- fread(paste0(resdir,
                            "svd.scalesnps.K100.F",fac, ".pseudoNTRUE.4400blocks.cell_type_results.txt"),
                     header=TRUE) %>%
      mutate(Factor = fac)
    svd_allres <- bind_rows(svd_allres, svdfile)
  }
  vb_allres %>% write_tsv(paste0(outdir, "vb.scaleTrue.K100.allfac.pseudoNTRUE.wt_cosine.4400blocks.cell_type_results.txt"))
  svd_allres %>% write_tsv(paste0(outdir,"svd.scalesnps.K100.allfac.pseudoNTRUE.4400blocks.cell_type_results.txt"))
}