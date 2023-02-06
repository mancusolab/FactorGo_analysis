#!/bin/bash
ldsc=/project/nmancuso_8/elezhang/software/ldsc/ldsc.py
cd /project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/variant/
P3dir="/project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3/plink_files"
annotdir="./Multi_tissue_gene_expr_1000Gv3_ldscores"
outdir="./baseline_v1.2"
snpdir="/project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/3_ldprune/sig_pval5e-08_nct2_maf1e-02_1kg/varinfo"

# pull the idx'th line from file params and store in variable named params
# params=`sed "${idx}q;d" params_chr`

# split strings in params variable by space and store in bash arg variables
# set -- junk $params
# shift

for chr in {1..22}
do
python ${ldsc} \
  --bfile ${P3dir}/1000G.EUR.QC.${chr} \
  --l2 \
  --ld-wind-cm 1.0 \
  --extract ${snpdir}/SNP1kg.${chr}.txt\
  --out ${outdir}/weights_factorgo.${chr}
done
