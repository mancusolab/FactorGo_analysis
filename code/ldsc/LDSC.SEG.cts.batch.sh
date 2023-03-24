#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=48:00:00
#SBATCH --mem-per-cpu=25Gb
#SBATCH --array=1-1025
#SBATCH -o ./Report/slurm.%a.SEG.vb.sep.out
#SBATCH --partition=main
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

# eg. start = 1, stop = 10
start=`python -c "print(1 + 20*int(int($idx-1)))"`
stop=$((start + 19))

ldsc=./software/ldsc/ldsc.py

cd ./data/PanUKBB/variant

P3dir="./DATA/ldsc/reference_files/1000G_EUR_Phase3/plink_files"
base_dir="./baseline_v1.2"
outdir="./results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/enrich/SEG/ldsc_cts"
weights_dir="./weights_factorgo"
sumstat_dir="./results/PanUKBB/sig_pval5e-08_nct2_maf1e-02_1kg/W_sumstats"

for IDX in `seq ${start} ${stop}`; do
  # run for FactorGo and tsvd separately for computational efficiency
  params=`sed "${IDX}q;d" param_SEG_file_annot_vb`
  echo "Running instance ${IDX} with params: ${params}"
  set -- junk $params
  shift

  sumstats=$1
  annot=$2

python ${ldsc} \
	--h2-cts ${sumstat_dir}/${sumstats}.gz\
	--w-ld-chr ${weights_dir}/weights_factorgo.\
	--ref-ld-chr ${base_dir}/baseline.\
	--overlap-annot\
        --n-blocks 4400 \
        --not-M-5-50 \
        --ref-ld-chr-cts ./SEG_ldcts/Annot${annot}.ldcts \
	--out ${outdir}/${sumstats}.Annot${annot}\
	--print-coefficients
done
