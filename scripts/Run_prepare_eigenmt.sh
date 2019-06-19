#!/bin/bash

## Script to prepare files for eigenMT
## Needs to be run once

## Usage: ./Run_prepare_eigenmt.sh

module load tabix

dir='../example'

chr='1 2 3 4 5 6 7 8 9 10 11 12 13 14 15 16 17 18 19 20 21 22'

tissue='GEUVADIS.445_samples'

# make directories
mkdir -p ${dir}/genotypes
mkdir -p ${dir}/positions
mkdir -p ${dir}/phenotypes

# 1) generate matrix eQTL format genotype files
vcf="${dir}/GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.vcf.gz"

for i in ${chr}
do
	echo "tabix -h ${vcf} chr${i} | python vcf2matrixeqtl.py | gzip -c > ${dir}/genotypes/chr${i}.txt.gz" | qsub -cwd -V -l mem=5G -N eigen_gt_${i}
	echo "tabix -h ${vcf} chr${i} | python vcf2matrixeqtl_pos.py | gzip -c > ${dir}/positions/chr${i}.txt.gz" | qsub -cwd -V -l mem=5G -N eigen_pos_${i}
done

# 2) generate phenotype position files - C position in bed files = beg + 1
for xtissue in ${tissue}
do
	echo "zcat ${dir}/${tissue}.expression.bed.gz | grep -v '#' | awk 'BEGIN{OFS=\"\t\";FS=\"\t\"}{if(NR==1){print(\"gene_id\",\"chrom_probe\",\"s1\",\"s2\")}; print(\$4,\$1,\$2+1,\$2+1)}' | gzip -c > ${dir}/phenotypes/${tissue}.txt.gz" | qsub -cwd -V -l mem=5G -N eigen_pheno
done
