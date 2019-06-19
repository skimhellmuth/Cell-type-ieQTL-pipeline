##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param gct
##' @param ...
##' @return
##' @author Sarah Kim-Hellmuth
inverse_quantile_normalization <- function(gct,
                                           ...) {
  gct = t(apply(gct, 1, rank, ties.method = "average", na.last="keep")); # rank row values, default of na.last is TRUE, where missing values in the data are put last
  gct = qnorm(gct / (ncol(gct)+1)); # with +1 inf values become 0, same for silva's approach
  return(gct)
}

##' ##' .. content for \description{} (no empty lines) ..
##'
##' .. content for \details{} ..
##' @title
##' @param dat
##' @param annot
##' @param tissue
##' @param celltype
##' @param out
##' @param outpath
##' @param verbose
##' @param ...
##' @return
##' @author Sarah Kim-Hellmuth
get.inv.norm <- function(dat,
                         annot=NULL,
                         tissue,
                         celltype=NULL,
                         out=F,
                         outpath=NULL,
                         verbose=T,
                         ...){
  sel <- which(annot$TISSUE_ABBRV==tissue)
  if (verbose) cat("No. of samples: ",length(sel),"\n")
  if (!is.null(annot)){
    id.s <- annot$SAMPID[sel]
    if (verbose) cat("No. of samples found: ",length(which(colnames(dat) %in% id.s)), "\n")
    scores.inv.norm <- inverse_quantile_normalization(dat[,id.s])
  } else {
    if (verbose) cat("Using all samples\n")
    scores.inv.norm <- inverse_quantile_normalization(dat)
  }
  if (!is.null(celltype)){
    scores.t <- as.data.frame(t(scores.inv.norm))
    scores.t$ID <- sub("^([^-]*-[^-]*).*", "\\1",rownames(scores.t))
    scores.inv.norm <- scores.t[,c("ID",celltype)]
  }
  scores.inv.norm <- scores.inv.norm[order(scores.inv.norm$ID),]
  if (out){
    tis.name <- unique(annot$TISSUE_NAME[sel])
    if(is.null(outpath)){
      outpath <- getwd()
    }
    if (!is.null(celltype)){
      if (verbose) cat("writing: ",paste0(outpath,"/", tis.name, "_normalized_xCell.", celltype,".txt \n"))
      write.table(scores.inv.norm, file=paste0(outpath,"/", tis.name, "_normalized_xCell.", celltype,".txt"), sep="\t",row.names=F, col.names=T, quote=F)
    } else {
      if (verbose) cat("writing: ",paste0(outpath,"/", tis.name, "_normalized_xCell.txt \n"))
      write.table(scores.inv.norm, file=paste0(outpath, "/", tis.name, "_normalized_xCell.txt"), sep="\t",row.names=T, col.names=T, quote=F)
    }
  } else {
    return(scores.inv.norm)
  }
}
