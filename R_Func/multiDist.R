#########################################################################################################################################
# This is function aims to compute multiple distance among sample from a matrix containing read counts                                                               #
# title: "'"                     #
# author:  "Quentin PETITJEAN[q.petitjean1@gmail.com]                                                                                   # 
# date: "15/06/2023"                                                                                                                    #
#########################################################################################################################################

multiDist <- function(dataList = NULL, # a matrix containing reads count, with sample name as column or a list of length 2, containing the read count and samples information
                      dist = NULL, # a vector of character string indicating the distances to compute and return
                      PhyloTree = NULL){ # an object of class phylo containing the phylogenetic tree correspoinding to the data or the path to the .nwk file containing the tree to import
  
  if (is.list(dataList)) {
    if (!"reads" %in% names(dataList)) {
      stop("dataList does not contain reads count, consider naming list element")
    }
  } else if (is.matrix(dataList)) {
    dataList <- list(reads = dataList)
  }
  
  # compute the specified distance
  Res <- list()
  vegdistList <- c("manhattan", "euclidean", "canberra", "clark", "bray", "kulczynski", "jaccard", "gower", "altGower", "morisita", "horn", "mountford", "raup", "binomial", "chao", "cao", "mahalanobis", "chisq", "chord", "hellinger", "aitchison", "robust.aitchison")
  for(i in dist){
    distance <- tolower(i)
    if(distance %in% vegdistList){
      Res[[i]] <- vegan::vegdist(t(dataList[["reads"]]), distance)
    }else if(distance %in% c("wunifrac",  "unifrac")){
      ## convert the data to Phyloseq object containing both OTU table and phylogenetic tree
      if(is.null(PhyloTree)){
        stop("a phylogenetic tree is needed to compute unifrac distances")
      } else if(class(PhyloTree) ==  "phylo"){
        PhyloTree <- PhyloTree
      } else if(is.character(PhyloTree)){
        PhyloTree <-
          ape::read.tree(PhyloTree)
      }
      Physeq <-
        phyloseq::phyloseq(phyloseq::otu_table(dataList[["reads"]], taxa_are_rows = T),
                           phyloseq::phy_tree(PhyloTree))
      Res[[i]] <- phyloseq::distance(Physeq, distance)
    } 
  }
  return(Res)
}
