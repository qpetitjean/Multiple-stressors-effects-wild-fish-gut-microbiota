#' Compute Multiple Beta Diversity Distances
#'
#' @description
#' The function \code{multiDist} calculates various beta diversity distance from a reads count dataset.
#' It supports standard distance metrics available in the \pkg{vegan} package (e.g., Bray-Curtis, Jaccard, Euclidean, etc.)
#' and also phylogenetic distances such as Unifrac and Weighted Unifrac, in case a phylogenetic tree is supplied.
#' The function returns a list of distance matrices, each corresponding to one of the specified metrics.
#'
#' @param dataList A numeric matrix containing reads count data with samples as columns, or a list 
#'   with a "reads" element containing such a matrix. If a matrix is provided, it is converted to a list.
#' @param dist A character vector specifying the distance metrics to compute. Supported metrics include standard
#'   distances such as "manhattan", "euclidean", "canberra", "bray", "jaccard", "hellinger", "aitchison", "robust.aitchison", etc.
#'   Additionally, "unifrac" and "wunifrac" are supported, which require a corresponding phylogenetic tree.
#' @param PhyloTree An optional parameter. A \code{phylo} object (from the \pkg{ape} package) representing the 
#'   phylogenetic tree corresponding to the taxa in the data, or a character string with the path to a Newick (.nwk) file.
#'   This parameter is required if "unifrac" or "wunifrac" distances are to be computed.
#'
#' @return A list of distance matrices. Each list element is named according to the corresponding distance metric.
#'
#' @details
#' The function first ensures that the input \code{dataList} is in list format (with a "reads" element). For standard
#' distances, it computes the distance on the transposed reads count matrix using \code{vegan::vegdist}. For phylogenetic
#' distances ("unifrac" or "wunifrac"), it requires a phylogenetic tree, which is used to create a \code{phyloseq} object,
#' and then computes the distance with \code{phyloseq::distance}.
#'
#' **Dependencies:**
#' This function depends on several packages. It internally calls functions from:
#' \itemize{
#'   \item \strong{vegan} (for computation of various distances),
#'   \item \strong{phyloseq} (for unifrac distances computation),
#'   \item \strong{ape} (for reading the phylogenetic tree)
#' }
#'
#' For reproducibility, please ensure that the required packages are installed and loaded.
#' Although not ideal to load packages within a function, the following libraries are explicitly 
#' loaded here to ensure that all dependent functions are available:
#'
#' @importFrom vegan vegdist
#' @importFrom phyloseq otu_table phy_tree phyloseq distance
#' @importFrom ape read.tree
#'
#' @author Quentin PETITJEAN [quentin.petitjean@inrae.fr]
#'
#' @date 15/06/2023
#'
#' @export

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
