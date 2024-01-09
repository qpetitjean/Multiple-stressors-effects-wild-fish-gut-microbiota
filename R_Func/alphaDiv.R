#############################################################################################################################
# This is function aims to compute alpha diversity indices from a matrix containing read counts per sample                  #
# Title: 
# author:  "Quentin PETITJEAN[q.petitjean1@gmail.com]                                                                       #
# date: "15/06/2023"                                                                                                        #                                                               #
#############################################################################################################################


alphaDiv <- function(Data = NULL, # a matrix containing reads count
                     indices = c("Richness", "Shannon", "Shannon.exp", "Simpson", "Chao1", "PD"), # a vector of character string indicating the alpha diversity indices to compute and return
                     PhyloTree = NULL # an object of class phylo containing the phylogenetic tree correspoinding to the data or the path to the .nwk file containing the tree to import
                     ){
  
  indices = tolower(indices)
  Res <- list()
  if ("richness" %in% indices) {
    # Compute Observed Richness
    Res[["Richness"]] <-
      apply(Data, 2, function(x)
        length(x[which(x != 0)]))
    Res[["Reads"]] <- apply(Data, 2, function(x)
      sum(x))
  }
  
  if ("shannon" %in% indices) {
    # Compute Shannon Index
    Res[["Shannon"]] <- vegan::diversity(t(Data), index = "shannon")
  }
  
  if ("shannon.exp" %in% indices) {
    if (!"shannon" %in% indices) {
      # Compute Shannon exponential
      Shannon <- vegan::diversity(t(Data), index = "shannon")
      Res[["Shannon.exp"]] <- exp(Shannon)
    } else{
      Res[["Shannon.exp"]] <- exp(Res[["Shannon"]])
    }
  }
  
  if ("simpson" %in% indices) {
    # Compute Simpson index
    Res[["Simpson"]] <- vegan::diversity(t(Data), index = "simpson")
  }
  
  if ("chao1" %in% indices) {
    # Compute Chao1
    if (!FALSE %in% sapply(Data, function(x)
      x == as.integer(x))) {
      Res[["Chao1"]] <- vegan::estimateR(t(Data))["S.chao1",]
      Res[["Chao1.se"]] <- vegan::estimateR(t(Data))["se.chao1",]
    } else{
      Res[["Chao1"]] <- rep(NA, ncol(Data))
      Res[["Chao1.se"]] <- rep(NA, ncol(Data))
    }
  }
  
  if ("pd" %in% indices) {
    if (is.null(PhyloTree)) {
      stop("a phylogenetic tree is needed to compute Faith' PD index")
    } else if (class(PhyloTree) ==  "phylo") {
      PhyloTree <- PhyloTree
    } else if (is.character(PhyloTree)) {
      PhyloTree <-
        ape::read.tree(PhyloTree)
    }
    # Compute Faith's PD
    Res[["PD"]] <-
      picante::pd(t(Data), PhyloTree, include.root = FALSE)[["PD"]]
  }
  
  # group all indices in a dataframe
  alphaDiv <- as.data.frame(do.call(cbind, Res))
  alphaDiv[["Samples"]] <- rownames(alphaDiv)
  rownames(alphaDiv) <- NULL
  
  return(alphaDiv)
}