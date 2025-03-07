# Script Title: Phylogenetic Tree Construction for 16s Sequences
#      
# Author: Quentin PETITJEAN
# Date Created: 03/2023
# Last Modified: 25/10/2023
# ==============================================================================
# Requirements: 
# - R version 4.2.3
# - Packages: 
#   - metabaR v1.0.0: For converting metabar data to FASTA format.
#   - Biostrings v2.66.0: For reading and writing DNA sequence data.
#   - DECIPHER v2.26.0: For aligning and adjusting DNA sequences.
#   - phangorn v2.11.1: For constructing and optimizing phylogenetic trees.
#   - ape v5.7-1: For tree manipulation and file output.
#   - ggtree v3.6.2: For phylogenetic tree visualization.
#   - ggplot2 v3.5.1: For enhanced plotting capabilities.
#   - viridis v0.6.3: For color scaling in plots.
# ==============================================================================
# Script Overview:
# This script builds a phylogenetic tree from 16s sequence data through the following steps:
# 1. Importing cleaned metabar data and converting it to FASTA format.
# 2. Aligning and adjusting the FASTA sequences.
# 3. Constructing an initial tree using the neighbor-joining method and evaluating multiple 
#    nucleotide substitution models.
# 5. Optimizing the maximum likelihood tree (including topology, branch lengths, and model parameters).
# 6. Assessing tree performance via AIC, bootstrap analyses, and various statistical tests.
# 7. Saving the final tree in Newick format and visualizing it with taxonomic annotations.
#
# ==============================================================================
# Usage:
# 1. Set the 'savingDir' variable to the directory containing your input files.
# 2. Ensure that the cleaned metabar data (fguts_Bact_agg_MergedRep.RDS) is available in the specified directory.
# 3. Run the script in an R environment to generate and visualize the phylogenetic tree.
# 4. The script outputs:
#    - A Newick tree file ("PhyloTree.nwk") containing the final tree.
#    - A visual tree plot ("phyloTree.svg") with taxonomic annotations.
# ==============================================================================

##############################################
#       	Install needed packages            #
##############################################

if(!require(metabaR)){
  install.packages("metabaR")
}
if(!require(Biostrings)){
  install.packages("Biostrings")
}
if(!require(DECIPHER)){
  install.packages("DECIPHER")
}
if(!require(phangorn)){
  install.packages("phangorn")
}
if(!require(ape)){
  install.packages("ape")
}
if(!require(ggtree)){
  install.packages("ggtree")
}
if(!require(ggplot2)){
  install.packages("ggplot2")
}
if(!require(viridis)){
  install.packages("viridis")
}

##########################################################
#  build the phylogenetic tree of 16s sequences          #
#########################################################

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported and saved
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs"

######################################################################################
# Convert the Metabar data to fasta sequences and aligned the sequences              #
######################################################################################

# import the cleaned data
metabarDat <- readRDS(file.path(savingDir, "Data/CleanedData", "fguts_Bact_agg_MergedRep.RDS"))

# save the dataset as fasta
metabaR::fasta_generator(metabarDat,
                         output_file = file.path(savingDir,"fguts_Bact_agg_Lab-Only.fasta"))

# import and align the dataset previously cleaned using the metabaR package and saved as .fasta
seqs <- Biostrings::readDNAStringSet(file.path(savingDir,"fguts_Bact_agg_Lab-Only.fasta"))

# Align the sequences 
AlSeqs <- DECIPHER::AlignSeqs(seqs)
AlSeqsAdj <- DECIPHER::AdjustAlignment(AlSeqs)
Biostrings::writeXStringSet(AlSeqsAdj, file.path(savingDir,"fguts_Bact_agg_Lab-Only_aligned_adj.fasta"))

# convert to PhyDat object
metabarDatFasta2 <- phangorn::read.phyDat(file.path(savingDir,"fguts_Bact_agg_Lab-Only_aligned_adj.fasta") ,format="fasta")

#########################
# Build the tree       #
########################

## uses DNA / AA sequences to compute pairwise distances 
distTree <- phangorn::dist.hamming(metabarDatFasta2)
treeIni <- phangorn::NJ(distTree)

# Compare different nucleotide or amino acid substitution models
ModComp <- phangorn::modelTest(
  tree = treeIni,
  metabarDatFasta2,
  model = "all",
  G = TRUE,
  I = TRUE,
  FREQ = TRUE
)

# retrieve the best one
bestMod <- phangorn::as.pml(ModComp,
                            rearrangement = "NNI")

# Now we can optimize the tree,
# using options to also optimize tree topology (optNni),
# base frequencies (optBf), and substitution rates (optQ).
# We will use a gamma distribution (optGamma) for variation in substitution
# rates at different sites in the sequence
OptiMod <-
  phangorn::optim.pml(
    bestMod,
    optNni = TRUE,
    optBf = TRUE,
    optQ = TRUE,
    optInv = TRUE,
    optGamma = TRUE,
    optEdge = TRUE,
    rearrangement = "NNI"
  )

# check whether optimized model is better than the raw one
## using Anova AIC and Shimodaira-Hasegawa Test
anova(bestMod, OptiMod)
phangorn::SH.test(bestMod, OptiMod)

## by checking tree distances and edge lenghts
phangorn::treedist(bestMod$tree, OptiMod$tree)
stats::setNames(data.frame(
  sum(treeIni$edge.length),
  sum(bestMod$tree$edge.length),
  sum(OptiMod$tree$edge.length)
),
c("TreeIni", "BestMod", "OptiMod"))


# Perform a bootstrap, starting with our 'best' maximum likelihood tree (OptiMod) with nearest neighbor interchanges
OptiModBoot <- phangorn::bootstrap.pml(OptiMod, bs=100, trees=TRUE, optNni=TRUE)

# we can overlay bootstrap supports on our maximum likelihood tree.
phangorn::plotBS(OptiMod$tree,OptiModBoot, type="phylogram",cex=0.5, bs.col = "firebrick", digits = 1)

# save the tree
PhyloTree <- OptiMod$tre
ape::write.tree(
  PhyloTree,
  file = file.path(savingDir,"PhyloTree.nwk"),
  append = FALSE,
  digits = 10,
  tree.names = FALSE
)

# visualize and save the tree (Supplementary material 2)
## retrieve taxonomic infos
OTUsNames <- data.frame(colnames(metabarDat$reads), metabarDat$motus$phylum, metabarDat$motus$family)

groupInfo <- split(PhyloTree$tip.label,
                   gsub("_\\w+", "", PhyloTree$tip.label))

## group family within phylum to display coloured path according to phylum
groups <- stats::setNames(lapply(unique(OTUsNames$metabarDat.motus.phylum), function(x)
  OTUsNames[grep(x, OTUsNames$metabarDat.motus.phylum), "metabarDat.motus.family"]
), unique(OTUsNames$metabarDat.motus.phylum))

## change the tip label (MOTUS_#) by family name
PhyloTree$tip.label <- OTUsNames[grep(paste(groupInfo$MOTU, collapse = "|"), OTUsNames$colnames.metabarDat.reads.),3]

## groups MOTUS according to the phylum 
PhyTree <- ggtree::groupOTU(PhyloTree, groups)

tree <- ggtree::ggtree(
  PhyTree,
  ggtree::aes(color = group),
  branch.length = 'none',
  layout = 'circular',
  ladderize = F,
  size = 1
) +
  ggtree::geom_tiplab(size = 2) +
  ggplot2::theme(
    legend.position = c(1, 1),
    legend.justification = c(0, 1),
    legend.title = ggplot2::element_text(size = 12),
    legend.text = ggplot2::element_text(size = 10)
  ) + 
  ggplot2::scale_colour_viridis_d(name="Phylum")
tree
ggplot2::ggsave(file=file.path(savingDir,"phyloTree.svg"), plot=tree, width=15, height=10,  dpi = 500)
