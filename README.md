# Cell type interaction-eQTL pipeline
#### A collection of scripts to estimate cell type estimates from bulk tissue gene expression data, perform interaction-eQTL (ieQTL) mapping and run colocalization analysis
---

## Dependencies

Following dependencies should be installed before running the pipeline:
- R (>= 3.2.2)
- Rpackages: xCell, GSEABase, coloc, [data.table](https://cran.r-project.org/package=data.table)
- Python (>= 3.6.4 and >=2.7.8)
- Python modules: pandas, numpy
- [tensorqtl](https://github.com/broadinstitute/tensorqtl) for eQTL mapping
- A modified [eigenMT](https://github.com/francois-a/eigenMT) python script to read tensorQTL output directly. It also include postprocessing steps that need [annotation.py](https://github.com/francois-a/rnaseq-utils) to produce the summary stat output format from GTEx v8. The original work can be found [here](https://www.sciencedirect.com/science/article/pii/S0002929715004929?via%3Dihub)
- tabix
- Example input data from the Geuvadis project are available [here](https://personal.broadinstitute.org/francois/geuvadis/)
---

## Workflow

To illustrate each step of this pipeline we will use data from the [Geuvadis](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/) project. Please download `GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.bed/.bim/.fam`, `GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.vcf.gz/.tbi`, `GEUVADIS.445_samples.expression.bed.gz/.tbi` and `GEUVADIS.445_samples.covariates.txt` from the link above and save in the `example` directory to run each of the scripts.

> **NOTE:** The purpose of this pipeline is to give an overview of analyses performed in the upcoming GTEx v8 cell type composition paper. It specifies individual parameters, provides helper functions and describes the environment each script was used in. To use it for your own data scripts will need to be adjusted accordingly. This workflow is not meant to be used as a standalone pipeline but as a template to reproduce GTEx v8 results or to use it in other projects.

#### 1. Estimate cell type abundance from bulk tissue gene expression data
The R-package [xCell](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1) calculates cell type enrichment scores from bulk tissue gene expression data. Following script `Run_xCell.sh` includes steps to read in the input file, add additional samples from GTEx v6p to ensure tissue-diversity, run xCell and outputs estimates of the cell type of interest to use as interaction covariate in tensorqtl.

```shell
./Run_xCell.sh
```
#### 2. Map cell type interaction-eQTLs

##### 2.1 Run eQTL interaction analysis
[Tensorqtl](https://www.biorxiv.org/content/10.1101/470138v3) can be used to map cis-eQTLs that interact with cell type estimates. Following script lists all tensorqtl parameters and parameters that are specific to the job scheduler slurm at NYGC.

```shell
celltype=B-cells; tissue=Cells_EBV-transformed_lymphocytes;
sbatch ./Run_tensorqtl_interaction.sh ${tissue} ${celltype}
```
##### 2.2 Correct for multiple testing
Permutationschemes adapted for interaction-eQTLs were not available when GTEx v8 cell type i-eQTLs were generated. We used eigenMT as alternative multiple testing correction method for cis-eQTL studies. Gene-level corrected p-values were further adjusted using Benjamini-Hochberg procedure.
> **NOTE:** To run `submit_eigenmt_py3.sh` genotype and phenotype input files need to be provided. These can be generated using the script `Run_prepare_eigenmt.sh`, which only needs to be run once.

```shell
celltype=B-cells; tissue=Cells_EBV-transformed_lymphocytes;
./Run_eigenmt.sh ${tissue} ${celltype}
```
The output file contains summary statistics for the interaction term, the genotype and cell type main effects, the eigenMT + Benjamini-Hochberg-corrected pvalues for the interaction term, and gene names and biotypes extracted from Gencode:
```
variant_id	gene_id	gene_name	biotype	phenotype_id	tss_distance	maf	ma_samples	ma_count	pval_g	b_g	b_g_se	pval_i	b_i	b_i_se	pval_gi	b_gi	b_gi_se	pval_emt	tests_emt	pval_adj_bh
chr1_940390_A_G_b38	ENSG00000227232.5	WASH7P	unprocessed_pseudogene	ENSG00000227232.5	910837	0.476119	475	638	0.497133	0.0285487	0.042019	0.0531376	0.252597	0.130363	9.70041e-05	-0.155152	0.0395342	0.0295863	305	0.259183
chr1_92750_G_A_b38	ENSG00000238009.6	RP11-34P13.7	lincRNA	ENSG00000238009.6	-36473	0.0701493	94	94	5.51215e-13	0.63642	0.0862951	0.18135	-0.175284	0.130988	4.50128e-05	0.385686	0.0938314	0.0178251	396	0.184471
```

#### 3. Perform colocalization analysis of cell type i-eQTLs and GWAS traits

##### 3.1 Extract significant cell type ieGenes
We ran colocalization analysis only on significant cell type ieGenes (FDR < 0.05). Following code extracts gene-ids of these ieGenes.
```shell
mkdir -p example/coloc
celltype=B-cells; tissue=Cells_EBV-transformed_lymphocytes;
cat example/${celltype}_ieqtl/eqtls_final/${tissue}_${celltype}_ieqtl.eigenMT.annotated.txt | awk '$21 < 0.05 {print $2}' > example/coloc/phenolist_${tissue}_${celltype}_ieqtl_fdr5.txt
```

Extract summary statistics of the entire cis-window of significant cell type ieGenes. This script needs about 60-80G of memory.
```shell
celltype=Neutrophils; tissue=Whole_Blood;
python3 extract_parquet2txt.py \
    -p ${celltype}_ieqtl/${tissue}*.parquet \
    -s phenolist_${tissue}_${celltype}_ieqtl_fdr5.txt \
    -o ${tissue}_${celltype}_ieqtl_fdr5.txt.gz \
    -d ${celltype}_ieqtl/subsets/
```

##### 3.1 Run coloc
```shell
celltype=Neutrophils; tissue=Whole_Blood;
./Run_coloc_gtex.v8i_114traits_prior_tq.sh ${tissue}_${celltype}_ieqtl_fdr5.txt.gz GTEx_Analysis_2017-06-05_v8_samplesize_abbrv_1.${celltype}.tissues.txt 1
```

##### 3.2 Postprocessing
```shell
celltype=Neutrophils; tissue=Whole_Blood;
./Run_coloc_post.sh coloc_${tissue}_${celltype}_ieqtl_fdr5_HEIGHT.Rda ${tissue}_${celltype}_ieqtl_fdr5 114traits
```

##### 3.3 Make locusplots
```shell
celltype=Neutrophils; tissue=Whole_Blood;
./Run_locusplot.sh Coloc.res4locuszoom_${tissue}_${celltype}_ieqtl_fdr5.txt 2 ieqtl locusplots/ieqtl
```
