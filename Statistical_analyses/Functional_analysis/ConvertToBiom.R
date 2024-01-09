##############################################################################################################
# Prepare file for conversion to BIOM format before using picrust2 for functional inference                  #
#############################################################################################################

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported 
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs"

# directory where the output files should be saved
DirSave <- file.path(savingDir, "BiomOutput")
if (length(list.dirs(DirSave)) == 0) {
  dir.create(DirSave)
}

#######################
# Import the data     #
#######################

# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "Varia_contam_fullv5.csv"), dec = ".", sep = ";")

# import the cleaned dataset (metabaR object) and select only gut samples
labFish <- readRDS(file.path(savingDir, "Preprocessing-Metabar/fguts_Bact_agg_MergedRep.RDS"))
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
source("https://raw.githubusercontent.com/itsmisterbrown/microfiltR/master/microfiltR_source_code.R")
library(phyloseq)

write.dataset(
  PhyloDat,
  paste0(DirSave, "/"),
  filePREFIX = "LabGm",
  writeFASTA = FALSE,
  rename = TRUE
)

# the next steps are performed using bash script in linux distribution