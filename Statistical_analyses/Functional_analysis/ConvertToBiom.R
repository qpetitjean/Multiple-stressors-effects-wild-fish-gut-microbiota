# Script Title: Prepare and Convert Metabar Object to BIOM Format for PICRUSt2 Functional Inference
#
# Author: Quentin PETITJEAN
# Date Created: 13/12/2023
# Last Modified: 13/12/2023
# ==============================================================================
# Requirements:
# - R version 4.2.3
# - Packages:
#   - metabaR v1.0.0: For processing and subsetting metabarcoding datasets.
#   - phyloseq v1.42.0: For constructing and manipulating phyloseq objects.
#   - seqinr v4.2-30: For writing FASTA files.
# - Source Files:
#   - microfiltR_source_code.R: Contains functions to convert phyloseq objects to BIOM format - from microfiltR (https://github.com/itsmisterbrown/microfiltR/tree/master): marker gene processing, filtering and exporting functions
# ==============================================================================
# Script Overview:
# This script prepares files for conversion to BIOM format, a necessary step before running
# PICRUSt2 for functional inference. The workflow includes:
# 1. Installing and loading required packages.
# 2. Specifying file paths for input data and output directories.
# 3. Importing experimental metadata and the cleaned metabaR object (filtered to include only fish gut samples).
# 4. Merging metadata with the metabaR object.
# 5. Converting the metabaR object into a phyloseq object (including OTU table, sample data, and taxonomy).
# 6. Generating a FASTA file containing sequence data, with sequence IDs prefixed by "ASV", formatted for PICRUSt2.
# 7. Converting the phyloseq object to a BIOM format file using a custom function from microfiltR.
# 8. Note: Subsequent steps for functional inference using PICRUSt2 are performed via a separate bash script.
# ==============================================================================
# Usage:
# 1. Ensure that all required source files and input data (e.g., DesignData.csv, fguts_Bact_agg_MergedRep.RDS)
#    are available in the specified directories.
# 2. Update the 'savingDir' variable to point to your data directory.
# 3. Run the script in an R environment to generate the FASTA and BIOM files.
# 4. Use the generated BIOM file with PICRUSt2 for functional inference.
# ==============================================================================

##############################################
#       	Install needed packages            #
##############################################
if(!require(metabaR)){
  install.packages("metabaR")
}
if(!require(phyloseq)){
  install.packages("phyloseq")
}
if(!require(seqinr)){
  install.packages("seqinr")
}

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported 
savingDir <- "W:/POSTDOC_INP_GOLFECH_2023/Outputs"

# directory where the output files should be saved
DirSave <- file.path(savingDir, "BiomOutput")
if (length(list.dirs(DirSave)) == 0) {
  dir.create(DirSave)
}

##############################################
#   Import some custom functions             #
##############################################
# the following function is imported from the following github repository: https://github.com/itsmisterbrown/microfiltR/tree/master
# and has been developed by Bryan Brown (itsmisterbrown github id)
source("https://raw.githubusercontent.com/itsmisterbrown/microfiltR/master/microfiltR_source_code.R") # import a function used to convert phyloseq object to Biom format

#######################
# Import the data     #
#######################

# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "DesignData/DesignData.csv"), dec = ".", sep = ";")

# import the cleaned dataset (metabaR object) and select only gut samples
labFish <- readRDS(file.path(savingDir, "Data/CleanedData", "fguts_Bact_agg_MergedRep.RDS"))
labFish <- metabaR::subset_metabarlist(labFish,
                                       table = "pcrs",
                                       indices = labFish$pcrs$matrix == "fish_gut")

# merge FullDat to the metabar sample object
labFish[["samples"]]["RowNames"] <- rownames(labFish[["samples"]])
labFish[["samples"]] <- merge(
  x = FullDat,
  y = labFish[["samples"]],
  by.x = "Ind",
  by.y = "Num_prlvt_Euth",
  all.y = TRUE
)
rownames(labFish[["samples"]]) <- labFish[["samples"]]$Ind
rownames(labFish[["reads"]]) <- labFish[["samples"]][["Ind"]][match(rownames(labFish[["reads"]]), labFish[["samples"]]$RowNames)]

# Convert to phyloseq object
## makes otu table for phyloseq
read.phylo <-
  phyloseq::otu_table(t(labFish[["reads"]]), taxa_are_rows = T)

## make sample table for phyloseq
sample.phylo <- phyloseq::sample_data(labFish[["samples"]])

##makes tax table for phyloseq
tax.phylo <-
  phyloseq::tax_table(as.matrix(labFish$motus[, c("phylum", "class", "order", "family", "genus", "species")]))

### merge phyloseq objects
PhyloDat  <- 
  phyloseq::phyloseq(read.phylo, sample.phylo, tax.phylo)

### Generate .fasta file containing sequences with Id starting by ">" (needed for picrust2)
seqinr::write.fasta(
  as.list(labFish$motus$sequence),
  paste0("ASV", seq_along(rownames(labFish$motus))),
  file.path(DirSave, "LabGmSeq.fasta"),
  as.string = TRUE
)

# convert sequences data to BIOM file (here we are not generating the FASTA using the write.dataset function because it does not returns sequences but instead OTUs Id)

library(phyloseq)
write.dataset(
  PhyloDat,
  paste0(DirSave, "/"),
  filePREFIX = "LabGm",
  writeFASTA = FALSE,
  rename = TRUE
)

# the next steps are performed using bash script in linux distribution (see Statistical_analyses/Functional_analysis/FuncInferences_Picrust2.txt)