#########################################################################################################################################
# This is function aims to compute multiple normalisation of reads count                                                                #
# title: "Multiple normalisation of 16s reads count data"
# author:  "Quentin PETITJEAN[q.petitjean1@gmail.com]                                                                                   # 
# date: "15/06/2023"                                                                                                                    #
#########################################################################################################################################

multiNorm <- function(dataList = NULL, # a matrix containing reads count, with sample name as column or a list of length 2, containing the read count and samples information
                      norm = c("original", "original.log",
                               "read.prop", "read.prop.log",
                               "read.CSS", "read.CSS.log",
                               "read.rare", "read.rare.log",
                               "read.deseq", "read.deseq.log",
                               "read.tmm", "read.tmm.log"), # a vector containing the list of normalization to compute and return, if a log transformation is specified the non-log version of the normalization should also be specified
                      TaxRow = F
                      ){
  if (is.list(dataList)) {
    if (!"reads" %in% names(dataList)) {
      stop("dataList does not contain reads count, consider naming list element")
    }
  } else if (is.matrix(dataList)) {
    dataList <- list(reads = dataList)
  }
  if(isFALSE(TaxRow)){
  dataList[["reads"]] <- t(dataList[["reads"]])
  }
  # compute the specified normalization
  Res <- list()
  # keep the original dataset and log2 transformation
  if ("original" %in% norm) {
    Res[["original"]] <- dataList[["reads"]]
  }
  if ("original.log" %in% norm) {
    Res[["original.log"]] <- log2(dataList[["reads"]] + 1)
  }
  
  ## Normalize the reads using various methods
  ### normalizes as proportions
  if ("read.prop" %in% norm) {
    read.prop <-
      stats::setNames(data.frame(matrix(
        ncol = ncol(dataList[["reads"]]),
        nrow = nrow(dataList[["reads"]])
      )), colnames(dataList[["reads"]]))
    for (i in seq(ncol(dataList[["reads"]]))) {
      read.prop[colnames(dataList[["reads"]])[i]] <-
        dataList[["reads"]][, i] / sum(dataList[["reads"]][, i]) * 10000
    }
    rownames(read.prop) <- rownames(dataList[["reads"]])
    Res[["read.prop"]] <- as.matrix(read.prop)
  }
  if ("read.prop.log" %in% norm) {
    Res[["read.prop.log"]] <- log2(Res[["read.prop"]] + 1)
  }
  
  ###normalizes as CSS
  if ("read.CSS" %in% norm) {
    perc <-
      metagenomeSeq::cumNormStatFast(obj = metagenomeSeq::newMRexperiment(dataList[["reads"]]),
                                     pFlag = F)
    Res[["read.CSS"]] <-
      metagenomeSeq::cumNormMat(metagenomeSeq::newMRexperiment(dataList[["reads"]]), perc, 1000)
  }
  
  if ("read.CSS.log" %in% norm) {
    Res[["read.CSS.log"]] <- log2(Res[["read.CSS"]] + 1)
  }
  
  #converts to phylo for rarefying and deseq
  if("read.rare" %in% norm | "read.rare.log" %in% norm | "read.deseq" %in% norm | "read.deseq.log" %in% norm){
    read.phylo <-
      phyloseq::otu_table(dataList[["reads"]], taxa_are_rows = T) #makes otu table for phyloseq
    #rownames(dataList[["samples"]]) <- dataList[["samples"]]$Num_prlvt_Euth
    sample.phylo <- phyloseq::sample_data(dataList[["samples"]])
    read.phylo  <- phyloseq::phyloseq(read.phylo, sample.phylo)
  }
  
  if ("read.rare" %in% norm | "read.rare.log" %in% norm) {
    ###normalizes by rarefying
    set.seed(2023)
    read.rare <- phyloseq::rarefy_even_depth(read.phylo)
    if("read.rare" %in% norm){
    Res[["read.rare"]] <- as(phyloseq::otu_table(read.rare), "matrix")
    }
  }
  if ("read.rare.log" %in% norm) {
    Res[["read.rare.log"]] <- log2(Res[["read.rare"]] + 1)
  }
  
  ###normalzes via deseq2
  if ("read.deseq" %in% norm | "read.deseq.log" %in% norm) {
    read.deseq <-
      phyloseq::phyloseq_to_deseq2(read.phylo, ~ pop) #converts to DESeq
    geoMeans <-
      apply(BiocGenerics::counts(read.deseq), 1, function(x, na.rm = TRUE) {
        exp(sum(log(x[x > 0]), na.rm = na.rm) / length(x))
      }) #calculates geometric means (taking into account the fact that there is zeros in the data)
    read.deseq <-
      DESeq2::estimateSizeFactors(read.deseq, geoMeans = geoMeans) #uses the means above to calaulte the size factors that are needed for variance stablization (this would not be neccissary if it were not for all the zeros)
    read.deseq <-
      DESeq2::varianceStabilizingTransformation(read.deseq, blind = F) #variance normalizes the data (this includes a log2(x+1) transformation)
    read.deseqTemp <-
      as.matrix(SummarizedExperiment::assay(read.deseq))
    if("read.deseq" %in% norm){
    Res[["read.deseq"]] <- (2 ^ read.deseqTemp) - 1
    Res[["read.deseq"]][Res[["read.deseq"]] < 0] <- 0
    }
  }
  if ("read.deseq.log" %in% norm) {
    Res[["read.deseq.log"]] <- read.deseqTemp
    Res[["read.deseq.log"]][Res[["read.deseq.log"]] < 0] <- 0
  }
  
  ###normalizes with TMM in edgeR
  if ("read.tmm" %in% norm) {
    tmm <- edgeR::calcNormFactors(dataList[["reads"]], method = "TMM")
    read.tmm <- NULL
    for (i in 1:length(tmm)) {
      col.i <- dataList[["reads"]][, i]
      col.i <- col.i / (tmm[i] * sum(col.i))
      read.tmm <- cbind(read.tmm, col.i)
    }
    colnames(read.tmm) <- colnames(dataList[["reads"]])
    Res[["read.tmm"]] <- read.tmm * 10000
  }
  
  if ("read.tmm.log" %in% norm) {
    Res[["read.tmm.log"]] <- log2(read.tmm + 1)
  }
  
  return(Res)
}
