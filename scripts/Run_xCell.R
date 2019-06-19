#!/usr/bin/Rscript

args <- commandArgs(TRUE);
INPUT <- args[1]
NAME <- args[2]
CORE <- args[3]
GTEX <- args[4] # e.g v6p_All_Tissues_gene_rpkm_4xCell_100samples.txt.gz
SAMPLEINFO <- args[5]
TISSUE <- args[6]
CELLTYPE <- args[7]

suppressPackageStartupMessages(library(xCell))
suppressPackageStartupMessages(library(GSEABase))
source("../Functions.R")

## read expression bed file, convert ensembl gene id to gene symbol and save file for future use
# library('biomaRt')
# mart <- useDataset("hsapiens_gene_ensembl", useMart("ENSEMBL_MART_ENSEMBL"))
# dat <- read.table(INPUT, head=T, comment.char="", as.is=T)
# dim(dat)
# dat$ens <- gsub("[.].*","",dat$gene_id)
# G_list <- getBM(filters= "ensembl_gene_id", attributes= c("ensembl_gene_id","hgnc_symbol"),values=dat$ens,mart= mart, uniqueRows=TRUE)
# dat.m <- merge(dat,G_list,by.x="ens",by.y="ensembl_gene_id", all=F, sort=F)
# dim(dat.m)
# dat.m <- dat.m[,c(2:4,ncol(dat.m),6:(ncol(dat.m)-1))] # bring back in bed format
# rem <- which(dat.m$hgnc_symbol=="")
# dat.m <- dat.m[-rem,]
# dat.m <- dat.m[!duplicated(dat.m$hgnc_symbol),]
# dim(dat.m)
# INPUT.NEW <- sub(".bed.gz","",INPUT)
# write.table(dat.m, file=gzfile(paste0(INPUT.NEW, "_hgnc.symbol.bed.gz")), col.names=T, row.names=F, quote=F, sep="\t") # can be uploaded in future runs

## read expression bed file that contains gene symbols instead of ensembl gene ids
## load 100 GTEx samples from v6p representing 41 tissues and add to existing matrix to ensure tissue diversity
dat <- read.table(INPUT, head=T, comment.char="", as.is=T)
dim(dat)
gtex <- read.table(GTEX, head=T, as.is=T)
dim(gtex)
dat.g <- merge(dat, gtex, by.x="hgnc_symbol", by.y="Description", all=F, sort=F)
dat.g <- dat.g[,c(2:4,1,5:ncol(dat.g))]
rownames(dat.g) <- dat.g$hgnc_symbol
dim(dat.g)

## run xCell
res <- xCellAnalysis(dat.g[,-c(1:4)], save.raw = T, file.name = NAME, parallel.sz=as.numeric(CORE)) # takes about 20min for v6p
dim(res)
res <- res[-c(65:67),] # delete cancer-related scores
saveRDS(res, file=paste0(NAME, ".rds"))

## get inverse normalized cell type estimates, reformat, and output as txt
annot <- read.table(SAMPLEINFO, head=T, as.is=T)
get.inv.norm(dat=res, annot=annot, tissue=TISSUE, celltype=CELLTYPE, out=T) # optional: add outputpath
