#!/bin/bash

## Script uses modified version of eigenmt to match francois' pipeline of generating final output files
## usage: module load python/3.7.1; trait=Neutrophils; for i in 1; do tissue=$(sed -n ''$i'p' tissues_${trait}.txt); qsub -cwd -V -j y -p -10 -N ei.${tissue} -pe smp 23 -l h_rt=18:00:00 -l mem=100G -m a -M skim@nygenome.org ./submit_eigenmt_py3.sh ${tissue} ${trait} ieqtl; done

set -eo pipefail

TISSUE=$1
CELLTYPE=$2
TYPE=$3
PREFIX=${TISSUE}_${CELLTYPE}_${TYPE}

MPATH=/gpfs/commons/groups/lappalainen_lab/skim/gtex_eqtl/v8/tensorqtl
FPATH=${CELLTYPE}_${TYPE}
EMTFPATH=/gpfs/commons/groups/lappalainen_lab/data/gtex/v8/eigenmt

## prepare directories
cd ${MPATH}/${FPATH}/

mkdir -p eqtls_final

python3 /gpfs/commons/groups/lappalainen_lab/skim/gtex_eqtl/v8/tensorqtl/run_eigenmt.py \
-q ${PREFIX}.cis_qtl_pairs.chr*.parquet \
-i /gpfs/commons/groups/lappalainen_lab/skim/gtex_eqtl/v8/tensorqtl/${CELLTYPE}_ieqtl/${TISSUE}*.eigenMT_input.txt.gz \
-g /gpfs/commons/groups/lappalainen_lab/data/gtex/v8/eigenmt/genotypes/* \
-gp /gpfs/commons/groups/lappalainen_lab/data/gtex/v8/eigenmt/positions/* \
-p /gpfs/commons/groups/lappalainen_lab/data/gtex/v8/eigenmt/phenotypes/${TISSUE}.txt.gz \
-o eqtls_final/ \
-x ${PREFIX} \
-s /gpfs/commons/groups/lappalainen_lab/skim/gtex_eqtl/v8/sample_lists/${TISSUE}.v8.samplelist.txt \
--parallel 23

exit
