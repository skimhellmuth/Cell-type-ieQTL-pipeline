#  -*- sh -*-
#$ -S /bin/bash
#$ -q res.q
#$ -cwd
#$ -V
#$ -l mem=12G
#$ -pe smp 4
#$ -j y
#$ -N  xcell

#set -eo
set -x

## usage: qsub ./Run_xCell.sh

module load R/3.2.2

input=GEUVADIS.445_samples.expression_hgnc.symbol.bed.gz
name=xCell_GEUVADIS.445_samples
core=4 # make sure it is the same as in qsub's -pe smp
gtex=v6p_All_Tissues_gene_rpkm_4xCell_100samples.txt.gz
sampleinfo=GEUVADIS.445_Sampleinfo.txt
tissue=LCL
celltype=B-cells

cd example

Rscript ../Run_xCell.R ${input} ${name} ${core} ${gtex} ${sampleinfo} ${tissue} ${celltype}
