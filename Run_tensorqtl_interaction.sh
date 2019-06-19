#!/bin/bash
#SBATCH --job-name=tq                  # Job name
#SBATCH --partition=gpu                      # Partition Name
#SBATCH --mail-type=FAIL                 # Mail events (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --mail-user=skim@nygenome.org        # Where to send mail
#SBATCH --mem=50gb                            # Job memory request
#SBATCH --time=4:00:00                       # Time limit 4 hours
#SBATCH --output=stdout_%j.log               # Standard output and error log
#SBATCH --gres=gpu:2       # To request for 1 GPU card

#set -eo
set -x

## log on to pe2 for gpu usage: sbatch ./Run_tensorqtl_interaction.sh Cells_EBV-transformed_lymphocytes B-cells

tissue=$1;
shift
celltype=$1;
shift

module unload R/3.4.1
module load cuda/10.0
source /gpfs/commons/groups/lappalainen_lab/software/anaconda3/bin/activate;
conda activate gpu;

p2f=example;
norm=normalized;
outdir="example/${celltype}_ieqtl"
prefix="${tissue}_${celltype}_ieqtl"

mkdir -p ${outdir}

/gpfs/commons/groups/lappalainen_lab/software/tensorqtl/tensorqtl/tensorqtl.py \
${p2f}/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered \
${p2f}/GEUVADIS.445_samples.expression.bed.gz \
${prefix} \
--covariates ${p2f}/GEUVADIS.445_samples.covariates.txt \
--mode cis_nominal \
--interaction ${p2f}/${tissue}_${norm}_xCell.${celltype}.txt \
--maf_threshold_interaction 0.05 \
-o ${outdir}
