library(tidyverse)
library(data.table)

args <- commandArgs(trailingOnly=TRUE)
datname <- args[1] # sig_pval5e-08_maf1e-02_nct2_1kg
snpfile <- args[2] # allvar.pruned
phefile <- args[3]d
options("scipen"=100, "digits"=4)

setwd(paste0("./data/PanUKBB/4_sumstats/", datname))

querylist_path <- paste0("./data/PanUKBB/3_ldprune/", datname, "/pruned/", snpfile)
manifest <- read_tsv(paste0("./data/PanUKBB/",phefile))
snplist <- read_tsv(querylist_path) %>% pull(rsid)

alldat <- data.table(rsid = snplist)
totalphe <- nrow(manifest)

for (idx in 1:totalphe){
  dat <- fread(paste0(manifest$filename[idx]))[, list(rsid, z_EUR)]
  names(dat)[names(dat) == "z_EUR"] <- paste0("z",idx)

  alldat <- alldat %>% left_join(dat, by = "rsid")

  if (idx %% 100 == 0){
    print(paste0("Finish on ", idx, "and current matrix size: ", dim(alldat)))
  }
}

## check colnames:
check_order <- sum(names(alldat)[2:(totalphe+1)] == paste0("z", 1:totalphe)) == totalphe
if (check_order == FALSE){
  alldat <- alldat %>% select(rsid, paste0("z", 1:totalphe))
}

alldat %>% fwrite(paste0("../", datname, "_allz.gz"), sep='\t')

### imputation with means
alldat <- alldat %>% as_tibble()
for (idx in 1:totalphe){
  onecol <- alldat[[idx+1]]
  if (sum(is.infinite(onecol)) > 0 | sum(is.na(onecol)) > 0){
    onecol[is.infinite(onecol)] <- NA
    onecol[is.na(onecol)] <- mean(onecol, na.rm=TRUE)
    alldat[[idx+1]] <- onecol
  }
  if (idx %% 500 == 0){
    print(paste0("Finish on impute", idx))
  }
}

check_NA_INF <- sum(alldat[, 2:(totalphe+1)])
if (is.na(check_NA_INF) == FALSE & is.infinite(check_NA_INF) == FALSE){
  alldat %>% fwrite(paste0("../", datname, "_allz.imputed.gz"), sep='\t')
}else{
  print("Data has NA or Inf")
}
