##############################################################################################################
#  Compute common taxonomic levels ratios indicating dysbiosis and test them                                 # 
# Bacteroidota/Proteobacteria ratio                                                                          #
# Firmicutes/Bacteroidota ratio                                                                              #
##############################################################################################################
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/DredgedLMM.R")  # import a function used to retrieve the summary, anova table and Rsquared of the most straightforward mixed effect model using rank selection (e.g., AICc)
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/multiNorm.R") # import a function used to compute multiple normalization methods
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GutMicrobiome/BxpltFunc.R")

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported 
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs"

###############################################################
# data preparation
###############################################################

# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "Varia_contam_fullv5.csv"), dec = ".", sep = ";")

# import the cleaned dataset (metabaR object) and select only gut samples
labFish <- readRDS(file.path(savingDir, "Preprocessing-Metabar/fguts_Bact_agg_MergedRep.RDS"))
labFish <- metabaR::subset_metabarlist(labFish,
                                       table = "pcrs",
                                       indices = labFish$pcrs$matrix == "fish_gut")

# retrieve OTUs taxonomy
OTUsTax <- data.frame(MOTU = colnames(labFish$reads),
                      phylum = labFish$motus$phylum, 
                      class = labFish$motus$class,
                      order = labFish$motus$order, 
                      family = labFish$motus$family, 
                      genus = labFish$motus$genus,
                      species = labFish$motus$species)

# aggregating the reads to the specified taxonomic level - family
taxLevel <- "phylum"
labFishPhy <- metabaR::aggregate_motus(labFish, groups = labFish$motus[[taxLevel]])

# merge FullDat to the metabar sample object
labFishPhy[["samples"]]["RowNames"] <- rownames(labFishPhy[["samples"]])
labFishPhy[["samples"]] <- merge(
  x = FullDat,
  y = labFishPhy[["samples"]],
  by.x = "Ind",
  by.y = "Num_prlvt_Euth",
  all.y = TRUE
)
rownames(labFishPhy[["samples"]]) <- labFishPhy[["samples"]]$RowNames


#######################################
# Normalize the data                  #
#######################################
normalizedDat <- multiNorm(
  dataList = labFishPhy,
  norm = c(
    "original.log"
  )
)

# Convert to phyloseq object 
## makes otu table for phyloseq
read.phylo <-
  phyloseq::otu_table(normalizedDat$original.log, taxa_are_rows = T)

## make sample table for phyloseq
sample.phylo <- phyloseq::sample_data(labFishPhy[["samples"]])

## merge phyloseq objects
PhyloDat  <- 
  phyloseq::phyloseq(read.phylo, sample.phylo)

# extract the count for selected phylums
Bacteroidota <- PhyloDat@otu_table[which(rownames(PhyloDat@otu_table) == "Bacteroidota")]
Firmicutes <- PhyloDat@otu_table[which(rownames(PhyloDat@otu_table) == "Firmicutes")]
Proteobacteria <- PhyloDat@otu_table[which(rownames(PhyloDat@otu_table) == "Proteobacteria")]

##########################################
# compute Firmicutes/Bacteroidota ratio  #
##########################################

FBr <- data.frame(Sample = colnames(Firmicutes), Firmicutes = as.vector(Firmicutes), Bacteroidota = as.vector(Bacteroidota))
FBr[["FBratio"]] <- FBr[["Firmicutes"]]/FBr[["Bacteroidota"]]
FBr <- cbind(FBr, as.data.frame(PhyloDat@sam_data))
FBr$FBratio[which(is.infinite(FBr$FBratio))] <- NA
FBr <- FBr[-which(is.na(FBr$FBratio)),]


# test the effect of treatment on the ratio
FBr[["Bac"]] <- as.character(FBr[["Bac"]])
Mod <- DredgedLMM(FBr, 
           expVar = "FBratio",
           respExpr = "(Contam + ContPop + Inj)^3 + Taille_mm_M",
           random = "(1 |Session / Bac) + (1 |Pop)",
           method = "weight",
           rank = "AICc")

performance::check_model(Mod$ModL)

newOrder <- c("NC", "C")
FBr[["ContamOrd"]] <- factor(FBr[["Contam"]], levels = newOrder)
Bxplt(FBr$ContamOrd, FBr$FBratio, 
      fill = c(
        adjustcolor("#006400", alpha.f = 0.3),
        adjustcolor("#B22222", alpha.f = 0.3)
      ),
      colpts = FBr[["Contam_col"]],
      boxwex = 0.5,
      ylab = "Firmicutes/Bacteroidota ratio",
      cex.axis = 1.1)

##############################################
# compute Bacteroidota/Proteobacteria ratio  #
#############################################
BPr <- data.frame(Sample = colnames(Proteobacteria), Proteobacteria = as.vector(Proteobacteria), Bacteroidota = as.vector(Bacteroidota))
BPr[["BPratio"]] <- BPr[["Bacteroidota"]]/BPr[["Proteobacteria"]]
BPr <- cbind(BPr, as.data.frame(PhyloDat@sam_data))
BPr$BPratio[which(is.infinite(BPr$BPratio))] <- NA
BPr <- BPr[-which(is.na(BPr$BPratio)),]

# test the effect of treatment on the ratio
BPr[["Bac"]] <- as.character(BPr[["Bac"]])
Mod <- DredgedLMM(BPr, 
                  expVar = "BPratio",
                  respExpr = "(Contam + ContPop + Inj)^3 + Taille_mm_M + Sexe",
                  random = "(1 |Session / Bac) + (1 |Pop)",
                  method = "weight",
                  rank = "AICc")

performance::check_model(Mod$ModL)

newOrder <- c("NC", "C")
BPr[["ContamOrd"]] <- factor(BPr[["Contam"]], levels = newOrder)
Bxplt(BPr$ContamOrd, BPr$BPratio, 
      fill = c(
        adjustcolor("#006400", alpha.f = 0.3),
        adjustcolor("#B22222", alpha.f = 0.3)
      ),
      colpts = BPr[["Contam_col"]],
      boxwex = 0.5,
      ylab = "Bacteroidota/Proteobacteria ratio",
      cex.axis = 1.1)
