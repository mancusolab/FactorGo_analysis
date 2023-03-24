## produce sumstats-like data for ldsc software:
## Header: SNP, CHR, BP, A1, A2, Z, N

library(tidyverse)
library(data.table)
options("scipen"=100, "digits"=4)

# datname="sig_pval5e-08_nct2_maf1e-02_1kg";K=100;method='vb';vbscale="True";svdscale='snps';pseudoN=TRUE;
args <- commandArgs(trailingOnly=TRUE)
datname <- args[1] # datname="sig_pval5e-08_nct2_maf1e-02_1kg"
K <- as.numeric(args[2])
method <- args[3] # svd or vb
vbscale <- args[4] # True or False
svdscale <- args[5] # snps, traits, none
pseudoN <- as.logical(args[6]) # whether to use Pseudo N or use posterior variance W

print(paste(datname, method, vbscale,svdscale, K, pseudoN))

setwd(paste0("/project/nmancuso_8/elezhang/projects/FA/results/PanUKBB/", datname))
outpath <- paste0("./W_sumstats/")

############# format loadings #############

# read in SNP chr, pos (ordered)
snpinfo <- read_tsv(paste0("/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/3_ldprune/",
                           datname, "/varinfo/allvar.pruned.1kgformat")) %>%
  select(SNP, CHR, BP, A1, A2)


# read Z and loading
if (method == "svd"){
  loading_path <- paste0("./results/PanUKBB/",
                         datname, "/tsvd.", datname, ".", svdscale,".k", K, ".V.tsv.gz")
  U_path <- paste0("./results/PanUKBB/",
                   datname, "/tsvd.", datname, ".", svdscale,".k", K, ".U.tsv.gz")
  D_path <- paste0("./results/PanUKBB/",
                   datname, "/tsvd.", datname, ".", svdscale,".k", K, ".D.tsv.gz")
  U_mat <- fread(U_path, header=FALSE)
  D_vec <- fread(D_path, header=FALSE) %>% pull(V1)
  Z2_mat <- (as.matrix(U_mat) %*% diag(D_vec))^2
  l2_mat <- fread(loading_path, header=FALSE) # pxk
  colnames(l2_mat) <- gsub("V", "X", colnames(l2_mat))
  print("read SVD loading and U matrix")
} else if (method == "vb"){
  loading_path <- paste0("./results/PanUKBB/",
                         datname, "/VB.", datname, ".k", K, ".scale", vbscale, ".Wm.tsv.gz")
  EZ2_path <- paste0("./results/PanUKBB/",
                    datname, "/VB.", datname, ".k", K, ".scale", vbscale, ".EZ2.tsv.gz")
  Zm_path <- paste0("./results/PanUKBB/",
                    datname, "/VB.", datname, ".k", K, ".scale", vbscale, ".Zm.tsv.gz")
  EZ2_mat <- fread(EZ2_path, header=FALSE)
  Zm_mat <- fread(Zm_path, header=FALSE) # nxk
  Z2_mat <- Zm_mat^2/(EZ2_mat - Zm_mat^2)
  l2_mat <- fread(loading_path, header=FALSE) # pxk
  colnames(l2_mat) <- gsub("V", "X", colnames(l2_mat))
  print(paste("read VB loading matrix and Z2 matrix from scaled", vbscale, "data"))
}else{
  print("svd or vb")
  quit()
}
print(dim(Z2_mat))
print(dim(l2_mat))

if (isTRUE(pseudoN)){
  sampleN <- read_tsv("./data/PanUKBB/4_sumstats/phe2483.SampleN.tsv", col_names = TRUE)

  # calculate cosine score (rowSums(cosine_phe) = 1), nxk
  cosine_phe  <- Z2_mat/rowSums(Z2_mat)

  # weight cosine score by sample size (kx1)
  Fn <- t(as.matrix(cosine_phe)) %*% as.matrix(sampleN)
}


####### Create sum stats file #######

create_file <- function(fac, Fn, W, snpinfo){

  sumstats <- snpinfo %>% mutate(Z = W[[fac]] * sqrt(Fn[[fac]]),
                                 N = Fn[[fac]])

  return(sumstats)
}

############# write out file for each line of param file #############

for (fac in 1:K){

  if (isTRUE(pseudoN)){
    dat <- create_file(fac, Fn, l2_mat, snpinfo)
  }

  outname <- ifelse(method == "vb",
                    paste0(method, ".scale", vbscale, ".K", K, ".F", fac, ".pseudoN", pseudoN, ".wt_cosine.gz"),
                    paste0(method, ".scale", svdscale, ".K", K, ".F", fac, ".pseudoN", pseudoN, ".gz"))

  if (sum(dat$Z < 0.001 * max(dat$N)) == nrow(dat)){
    dat %>% write_tsv(paste0(outpath, outname))
    print(paste0("Finished: ", fac))
  }else{
    print("has outlier Z that will be removed from ldsc")
  }
}
