### Count number sig var across files in
### TODO: generate list of variants other than this filtering, so that we can count for those variants
library(tidyverse)
library(data.table)

logpval_threshold <- log(5e-8)
outpath <- "./variant/"

indir <- "./PanUKBB/"
setwd(indir)

manifest_list <- read_tsv("./phe2483.PanUKBB.tsv") %>% pull(filename)

## Read in each sig file
## 1. get variant list from each phe file
vardf <- tibble()
for (idx in seq_along(manifest_list)){
  onefile <- read_tsv(paste0("./1_sigSNP/", manifest_list[[idx]]), col_types = 'cdccd') %>%
    filter(pval_EUR < logpval_threshold) %>%
    select(chr:alt)

    vardf <- bind_rows(vardf, onefile)
    print(idx)
}
vardf <- vardf %>% distinct(chr,pos,ref,alt)

# read SNP from SNP list that has 5e-8 signal
# vardf <- read_tsv("./variant/SNPlist_hassig_5e-8.tsv.gz", col_types = 'cdccd')
# snplist <- read_tsv("./3_ldprune/sig_pval5e-08_nct2_maf1e-02_1kg/varinfo/allvar.pruned", col_types = 'cdccccd')

# count for filtered variants: INFO>0.9, maf>0.01,high quality, not low confidence
snplist <- read_tsv("./2_varqc/full_variant_info0.9_maf1e-02_highqualTRUE.gz", col_types = 'cdccccd')
vardf <- snplist %>% left_join(vardf, by=c("chr", "pos", "ref", "alt"))

### Count
ctsig_trait_df <- tibble(filename = manifest_list, ct_sig = 0)
ctsig_var_df <- vardf %>% mutate(ct_sig = 0)

for (idx in seq_along(manifest_list)){
  onefile <- read_tsv(paste0("./1_sigSNP/", manifest_list[[idx]]), col_types = 'cdccd') %>%
    filter(pval_EUR < logpval_threshold)

  # count for variants INFO>0.9, maf>0.01,high quality, not low confidence & pval < _
  onefile_joined <- vardf %>% left_join(onefile, by=c("chr", "pos", "ref", "alt")) %>%
    mutate(is_sig = ifelse(is.na(pval_EUR), 0, 1))

  ctsig_trait_df$ct_sig[[idx]] <- sum(onefile_joined$is_sig)
  ctsig_var_df$ct_sig <- ctsig_var_df$ct_sig + onefile_joined$is_sig
  print(idx)
}

## Count number of sig traits per variant
ctsig_trait_df %>% fwrite("./variant/Trait_count_sig_5e-8_filtered.tsv.gz", sep="\t")

# Count number of sig var per trait
ctsig_var_df %>% fwrite("./variant/Var_count_sig_5e-8_filtered.tsv.gz", sep="\t")
