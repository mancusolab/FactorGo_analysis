#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=24:00:00
#SBATCH --mem-per-cpu=25Gb
#SBATCH --array=1-200
#SBATCH -o ./Report/slurm.%a.randomgene.sep.out
#SBATCH --partition=main
#SBATCH --account=nmancuso_8
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

ldsc=/project/nmancuso_8/elezhang/software/ldsc/ldsc.py

cd /project/nmancuso_8/elezhang/projects/FA/data/PanUKBB/variant

P3dir="/project/gazal_569/DATA/ldsc/reference_files/1000G_EUR_Phase3/plink_files"
base_dir="./baseline_v1.2"
outdir="/project/nmancuso_8/elezhang/projects/FA/results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/enrich/random_gene_annotations/ldsc_cts"
weights_dir="./weights_factorgo"
sumstat_dir="/project/nmancuso_8/elezhang/projects/FA/results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/W_sumstats"

# pull the idx'th line from file params and store in variable named params
params=`sed "${idx}q;d" param_pLI10_file_annot`

# split strings in params variable by space and store in bash arg variables
set -- junk $params
shift
  
sumstats=$1

python ${ldsc} \
	--h2-cts ${sumstat_dir}/${sumstats}.gz\
	--w-ld-chr ${weights_dir}/weights_factorgo.\
	--ref-ld-chr ${base_dir}/baseline.\
	--overlap-annot\
        --n-blocks 4400 \
        --not-M-5-50 \
        --ref-ld-chr-cts ./random_gene_annotations/random_gene_annotations.ldcts \
	--out ${outdir}/${sumstats}\
	--print-coefficients
