#' Compute Alpha Diversity Indices from Reads Count Data
#'
#' @description
#' The function \code{alphaDiv} computes multiple alpha diversity indices from a matrix of reads counts.
#' It calculates observed richness, Shannon diversity (and its exponential form), Simpson index,
#' Chao1 richness estimator (with standard error), and Faith's Phylogenetic Diversity (PD) in case a phylogenetic tree is supplied.
#' The results are returned as a data frame with one row per sample.
#'
#' @param Data A numeric matrix containing reads count per sample. Taxa are represented by rows and samples by columns.
#' @param indices A character vector indicating which alpha diversity indices to compute.
#'        Default is \code{c("Richness", "Shannon", "Shannon.exp", "Simpson", "Chao1", "PD")}. 
#'        Index names are case-insensitive.
#' @param PhyloTree An optional parameter. A \code{phylo} object (from the \pkg{ape} package) representing the phylogenetic tree
#'        corresponding to the taxa in \code{Data}, or a character string providing the path to a Newick (.nwk) file.
#'        This parameter is required if "PD" (Faith's Phylogenetic Diversity) is specified in `indices`
#'
#' @return A data frame containing the computed alpha diversity indices for each sample. The output includes:
#' \describe{
#'   \item{Richness}{Observed richness: the number of taxa with non-zero counts per sample.}
#'   \item{Shannon}{Shannon diversity index.}
#'   \item{Shannon.exp}{The exponential of the Shannon index.}
#'   \item{Simpson}{Simpson diversity index.}
#'   \item{Chao1}{The Chao1 richness estimator.}
#'   \item{Chao1.se}{The standard error of the Chao1 estimator.}
#'   \item{PD}{Faith's Phylogenetic Diversity (if a phylogenetic tree is provided).}
#'   \item{Reads}{The total reads per sample.}
#'   \item{Samples}{Sample identifiers.}
#' }
#'
#' @details
#' The function transforms the input \code{Data} matrix as necessary and uses functions from the \pkg{vegan} package
#' to compute Shannon and Simpson indices as well as the Chao1 estimator.
#' If "PD" is requested in \code{indices}, a phylogenetic tree must be provided via the \code{PhyloTree} parameter.
#' The function accepts the tree either as a \code{phylo} object or as a file path to a Newick formatted tree,
#' which is then imported using \code{ape::read.tree}.
#'
#' **Dependencies:**
#' This function depends on several packages. It internally calls functions from:
#' \itemize{
#'   \item \strong{vegan} (for computation of various indices),
#'   \item \strong{picante} (for faith phylogenetic distance computation),
#'   \item \strong{ape} (for reading the phylogenetic tree)
#' }
#'
#' For reproducibility, please ensure that the required packages are installed and loaded.
#' Although not ideal to load packages within a function, the following libraries are explicitly 
#' loaded here to ensure that all dependent functions are available:
#'
#' @importFrom vegan diversity estimateR
#' @importFrom picante pd
#' @importFrom ape read.tree
#'
#' @author Quentin PETITJEAN [quentin.petitjean@inrae.fr]
#'
#' @date 15/06/2023
#'
#' @export

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