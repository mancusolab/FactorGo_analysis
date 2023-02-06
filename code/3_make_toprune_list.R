###### prepare Hail LD prune input ######
## Remove indels, multi-allelic
## Keep variants only with LDSC annotations (1000GP3)

library(tidyverse)
library(data.table)
options("scipen"=100, "digits"=4)

args <- commandArgs(trailingOnly=TRUE)
datname <- args[1] # datname="sig_pval5e-08_nct2_maf1e-02";outdir_name="sig_pval5e-08_nct2_maf1e-02_1kg"
outdir_name <- args[2] # outdir_name="sig_pval5e-08_nct2_maf1e-02_1kg"

indir <- "/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/2_varqc/"
outdir <- paste0("/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/3_ldprune/", outdir_name, "/toprune/")
# outdir <- paste0("/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/3_ldprune/", outdir_name, "/varinfo/hg38/")
P3dir <- "/project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3/plink_files/" # 1000G.EUR.QC.{chr}.bim
# plink_hg38 <- "/project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3_hg38/plink_files/"

allvar <- fread(paste0(indir, datname, ".gz"))

# check allele flip
allele.qc <- function(a1,a2,ref1,ref2) {
  # a1 and a2 are the first data-set
  # ref1 and ref2 are the 2nd data-set
  # Make all the alleles into upper-case, as A,T,C,G:
  a1 = toupper(a1)
  a2 = toupper(a2)
  ref1 = toupper(ref1)
  ref2 = toupper(ref2)
  # Strand flip, to change the allele representation in the 2nd data-set
  strand_flip = function(ref) {
    flip = ref
    flip[ref == "A"] = "T"
    flip[ref == "T"] = "A"
    flip[ref == "G"] = "C"
    flip[ref == "C"] = "G"
    flip
  }
  flip1 = strand_flip(ref1)
  flip2 = strand_flip(ref2)
  snp = list()
  # not doing this: Remove strand ambiguous SNPs (scenario 3)
  #snp[["keep"]] = !((a1=="A" & a2=="T") | (a1=="T" & a2=="A") | (a1=="C" & a2=="G") | (a1=="G" & a2=="C"))
  snp[["keep"]] = rep(TRUE, length(a1))
  # Remove non-ATCG coding
  snp[["keep"]][ a1 != "A" & a1 != "T" & a1 != "G" & a1 != "C" ] = F
  snp[["keep"]][ a2 != "A" & a2 != "T" & a2 != "G" & a2 != "C" ] = F
  # as long as scenario 1 is involved, sign_flip will return TRUE
  snp[["sign_flip"]] = (a1 == ref2 & a2 == ref1) | (a1 == flip2 & a2 == flip1)
  # as long as scenario 2 is involved, strand_flip will return TRUE
  snp[["strand_flip"]] = (a1 == flip1 & a2 == flip2) | (a1 == flip2 & a2 == flip1)
  # remove other cases, eg, tri-allelic, one dataset is A C, the other is A G, for example.
  exact_match = (a1 == ref1 & a2 == ref2) 
  snp[["keep"]][!(exact_match | snp[["sign_flip"]] | snp[["strand_flip"]])] = F
  return(snp)
}

## read in 1kg variants
for (idx in 1:22){
  varfile <- allvar[chr == idx]
  if (idx == 6){
    ## exclude MHC regions
    varfile <- varfile[(pos < 25000000 | pos > 34000000)]
  }
  
  P3file <- fread(paste0(P3dir, "1000G.EUR.QC.", idx, ".bim"))
  # P3file <- fread(paste0(plink_hg38, "1000G.EUR.hg38.", idx, ".bim"))
  names(P3file) <- c('chr', 'SNP', 'cM', 'pos', 'A1', 'A2')
  to_check <- varfile %>% inner_join(P3file, by = c("chr", "pos"))
  
  P3file_checkindel <- P3file %>% filter((str_length(A1) > 1 | str_length(A2) > 1))
  if (nrow(P3file_checkindel) < 1){
    to_keep <- allele.qc(to_check$ref, to_check$alt, to_check$A1, to_check$A2)
    print("no indels")
  }else{
    break
  }
  
  ## allele matching
  outvar <- to_check[to_keep$keep]
  
  # for splitted variants, remove multi-allelic variants
  outvar %>%
    group_by(chr, pos, ref) %>%
    add_count() %>%
    filter(n == 1) %>%
    select(-n) %>%
    ungroup() %>%
    unite("variant", chr:alt, sep=":") %>%
    select(variant, maf_EUR) %>%
    fwrite(paste0(outdir, "chr", idx, ".toprune"), sep="\t")

  check_indel <- outvar %>% filter((str_length(ref) > 1 | str_length(alt) > 1))
  check_snpid <- length(outvar$rsid) == length(unique(outvar$rsid))
  print(check_indel)
  print(check_snpid)
}
