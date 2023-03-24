#!/bin/bash
#SBATCH --ntasks=1
#SBATCH --time=1:00:00
#SBATCH --mem-per-cpu=55Gb
#SBATCH --gres=gpu:2
#SBATCH --array=1
#SBATCH --mail-type=all
#SBATCH --mail-user=zzhang39@usc.edu

if [ ! $SLURM_ARRAY_TASK_ID ]; then
    idx=$1
else
    idx=$SLURM_ARRAY_TASK_ID
fi

cd  .results/PanUKBB/
dat="./data/PanUKBB/4_sumstats"

# configure
# module load gcc
module load gcc/8.3.0
export CC=`which gcc`
export OBJCC=`which gcc`
export CXX=`which g++`
export CPP=`which cpp`
export GCC=`which gcc`
export AR=`which ar`
export RANLIB=`which ranlib`
export NM=`which nm`
export GCC_AR=`which ar`
export GCC_RANLIB=`which ranlib`
export GCC_NM=`which nm`

module load cuda/11.2.0
export XLA_FLAGS=--xla_gpu_cuda_data_dir=/spack/apps/linux-centos7-x86_64/gcc-8.3.0/cuda-11.2.0-7uwimxj27s4cptafhkw6a6fpyqf5nw4c/
export LD_LIBRARY_PATH=/spack/apps/linux-centos7-x86_64/gcc-8.3.0/cuda-11.2.0-7uwimxj27s4cptafhkw6a6fpyqf5nw4c/lib64/:$LD_LIBRARY_PATH

# pull the idx'th line from file params and store in variable named params
params=`sed "${idx}q;d" params_vb`

# split strings in params variable by space and store in bash arg variables
set -- junk $params
shift

elbo_tol=0.001
max_run=20000000

datfile=$1
k=$2
scale=$3
platform=$4
Nfile="phe2483.SampleN.tsv"
seed=1

python runVB.py \
    ${dat}/${datfile}_allz.imputed.gz \
    ${dat}/${Nfile} \
    -k $k \
    --hyper 1e-5 1e-5 1e-5 1e-5 1e-5\
    --rate 1\
    --scaledat $scale \
    --max-iter $max_run \
    -o ./${datfile}/VB.${datfile}.k${k}.scale${scale}\
    -p $platform
