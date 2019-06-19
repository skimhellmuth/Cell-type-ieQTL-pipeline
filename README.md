# Cell type interaction eQTL pipeline
#### A pipeline to estimate cell type estimates from bulk tissue gene expression data, perform ieQTL mapping and run colocalization analysis
---

## Dependencies

Following dependencies should be installed before running the pipeline:
- Rpackages: xCell, GSEABase, coloc, data.table
- python
- [tensorqtl](https://github.com/broadinstitute/tensorqtl) for eQTL mapping
- A modified [eigenMT](https://github.com/francois-a/eigenMT) python script to read tensorQTL output directly. It also include postprocessing steps that need [annotation.py](https://github.com/francois-a/rnaseq-utils) to produce the summary stat output format from GTEx v8. The original work can be found [here](https://www.sciencedirect.com/science/article/pii/S0002929715004929?via%3Dihub)
- Example input data from the Geuvadis project are available [here](https://personal.broadinstitute.org/francois/geuvadis/)
---

## Basic usage

Please download `GEUVADIS.445_samples.GRCh38.20170504.maf01.filtered.bed/.bim/.fam`, `GEUVADIS.445_samples.expression.bed.gz/.tbi` and `GEUVADIS.445_samples.covariates.txt` from the link above and save in the example directory to be able to run the pipeline. Scripts need to be adjusted to use with your own data.

> **NOTE:** Data from the [Geuvadis](https://www.ebi.ac.uk/arrayexpress/experiments/E-GEUV-1/) project will be used for illustration purposes of the workflow. Results (e.g. cell type estimates generated from LCLs) are not meant to be biologically interpreted.

#### 1. Estimate cell type abundance from bulk tissue gene expression data
The R-package [xCell](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-017-1349-1) calculates cell type enrichment scores from bulk tissue gene expression data. Following script `Run_xCell.sh` includes steps to read in the input file, add additional samples from GTEx v6p to ensure tissue-diversity, run xCell and outputs estimates of the cell type of interest to use as interaction covariate in tensorqtl.

```shell
./Run_xCell.sh
```
#### 2. Map cell type interaction-eQTLs

##### 2.1 Run eQTL interaction analysis
To run `Run_tensorqtl_interaction.sh` tensorqtl

```shell
sbatch ./Run_tensorqtl_interaction.sh Cells_EBV-transformed_lymphocytes B-cells
```
##### 2.2 Correct for multiple testing

```shell
./submit_eigenmt_py3.sh Cells_EBV-transformed_lymphocytes B-cells
```

##### 2.3 Extract significant hits

#### 3. Perform colocalization analysis of cell type i-eQTLs and GWAS traits
