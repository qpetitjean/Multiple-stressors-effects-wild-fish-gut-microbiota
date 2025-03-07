# Script Title: Compute Taxonomic Ratios as Dysbiosis Indicators and Test Treatment Effects in Fish Gut Microbiota
#
# Author: Quentin PETITJEAN
# Date Created: 06/05/2024
# Last Modified: 06/05/2024
# ==============================================================================
# Requirements:
# - R version 4.2.3
# - Packages:
#   - metabaR v1.0.0: For processing, subsetting, and aggregating metabarcoding datasets.
#   - phyloseq v1.42.0: For constructing and manipulating phyloseq objects.
#   - performance v0.10.3: For checking model performance.
#   - kableExtra v1.3.4: For creating enhanced HTML summary tables.
#
# - Source Files:
#   - multiNorm.R: Function to compute multiple normalization methods.
#   - DredgedLMM.R: Function to perform model dredging for LMM and retrieve summaries.
#   - BxpltFunc.R: Function to generate customized boxplots.
# ==============================================================================
# Script Overview:
# This script computes two taxonomic ratios that are used as indicators of dysbiosis:
#  1. The Firmicutes/Bacteroidota ratio
#  2. The Bacteroidota/Proteobacteria ratio
#
# For each ratio, the workflow is as follows:
# 1. Import the experimental metadata and the cleaned metabaR object (focusing on gut samples).
# 2. Retrieve OTU taxonomy and aggregate the read counts at the phylum level.
# 3. Normalize the aggregated data using multiple normalization methods (via the multiNorm function).
# 4. Extract the read counts for the phyla of interest.
# 5. Compute the ratios (Firmicutes/Bacteroidota and Bacteroidota/Proteobacteria) for each sample.
# 6. Test the effect of treatments (and covariates) on these ratios using a linear mixed model 
#    framework with model dredging (via the DredgedLMM function).
# 7. Generate diagnostic plots for model performance and significant effects; save these plots as TIF files.
# 8. Save the statistical test outputs as RDS files for later reference.
# 9. Compiling and exporting HTML summary tables of statistical results for each taxonomic ratios.
#
# Usage:
# 1. Ensure that all required source files and input data files (e.g., DesignData.csv, fguts_Bact_agg_MergedRep.RDS)
#    are available in the specified directories.
# 2. Update the 'savingDir' variable to the correct directory path for your data and outputs.
# 3. Run the script in an R environment to compute the taxonomic ratios, test treatment effects, and generate 
#    diagnostic plots and summary tables.
# 4. The outputs include:
#    - TIF files of model diagnostic plots for each ratio under different normalization methods.
#    - RDS files containing the LMM results for the taxonomic ratios.
#    - HTML tables summarizing statistical results
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
if(!require(performance)){
  install.packages("performance")
}
if(!require(kableExtra)){
  install.packages("kableExtra")
}

#######################
# specify some path   #
#######################
# the path to the directory from where file should be imported 
savingDir <- "W:/POSTDOC_INP_GOLFECH_2023/Outputs"
DirSave <- file.path(savingDir, "Normalized_Data")

##############################################
#   Import some custom functions             #
##############################################
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/multiNorm.R") # import a function used to compute multiple normalization methods
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/DredgedLMM.R")  # import a function used to retrieve the summary, anova table and Rsquared of the most straightforward mixed effect model using rank selection (e.g., AICc)
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/BxpltFunc.R") # import a function compute draw customized boxplots

#######################
# Import the data     #
#######################

# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "Data/DesignData", "DesignData.csv"), dec = ".", sep = ";")

# import the cleaned dataset (metabaR object) and select only gut samples
labFish <- readRDS(file.path(savingDir, "Data/CleanedData", "fguts_Bact_agg_MergedRep.RDS"))
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
    "original",
    "original.log",
    "read.prop",
    "read.prop.log",
    "read.CSS",
    "read.CSS.log",
    #"read.rare",
    #"read.rare.log",
    "read.deseq",
    "read.deseq.log",
    "read.tmm",
    "read.tmm.log"
  )
)

#####################################################################
# Compute and test treatement effects on the ratios                #
####################################################################

RatioRes <- list()

for (i in names(normalizedDat)) {
  
  # Convert to phyloseq object 
  ## makes otu table for phyloseq
  read.phylo <-
    phyloseq::otu_table(normalizedDat[[i]], taxa_are_rows = T)
  
  ## make sample table for phyloseq
  sample.phylo <- phyloseq::sample_data(labFishPhy[["samples"]])
  
  ## merge phyloseq objects
  PhyloDat  <- 
    phyloseq::phyloseq(read.phylo, sample.phylo)
  
  # extract the count for selected phylums
  Bacteroidota <- PhyloDat@otu_table[which(rownames(PhyloDat@otu_table) == "Bacteroidota")]
  Firmicutes <- PhyloDat@otu_table[which(rownames(PhyloDat@otu_table) == "Firmicutes")]
  
 
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
  BestMod <- DredgedLMM(FBr, 
                    expVar = "FBratio",
                    respExpr = "(Contam + ContPop + Inj)^3 + Taille_mm_M",
                    random = "(1 |Session / Bac) + (1 |Pop)",
                    method = "weight",
                    rank = "AICc")

  if (length(list.dirs(file.path(DirSave, "Signif_Effects"))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects"))
  }
  if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio"))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio"))
  }
  if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio", i))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio", i))
  }
  if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio", i, "FBratio"))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio", i, "FBratio"))
  }
  toSave <- file.path(DirSave, "Signif_Effects", "TaxRatio", i, "FBratio", paste0("Mod_Perf_", "FBratio","_",  i, ".tif" ))
  if(file.exists(toSave)){
    unlink(toSave)
  } 
  
  grDevices::tiff(filename = toSave,
                  width = 1200,
                  height = 1200,
                  res = 100,
                  compression = "lzw",
                  pointsize = 8)
  print(performance::check_model(BestMod[["ModL"]]))
  dev.off()
  
  # store the result of the model
  RatioRes[[i]][["FBRatio"]] <- BestMod
  
  # display the significant results 
  toKeep <- rownames(BestMod$ModAnov)[-1]
  
  if(length(toKeep) == 0){
    next
  }
  ##### save the figs representing difference of PCOA Scores (composition among significant treatments)
  toSave <- file.path(DirSave, "Signif_Effects", "TaxRatio", i, "FBratio", paste0("SignifPlot", "FBratio","_",  i, ".tif" ))
  if(file.exists(toSave)){
    unlink(toSave)
  } 
  grDevices::tiff(filename = toSave,
                  width = 1200,
                  height = 1200,
                  res = 300,
                  compression = "lzw",
                  pointsize = 8)
  par(mfrow = c(ifelse(ceiling(sqrt(length(toKeep)))^2 > length(toKeep), ceiling(sqrt(length(toKeep))),2), ceiling(sqrt(length(toKeep)))), 
      mar = c(4, 4, 4, 2))
  for(j in seq_along(toKeep)){ 
    if(!length(grep(":", toKeep[[j]])) > 0){
      keepVar <- toKeep[[j]]
    }else{
      keepVar <-  gsub(":", ".", toKeep[[j]])
    }
    FBrOrd <- FBr[order(FBr[[keepVar]]), ]
    if(is.character(FBrOrd[[keepVar]]) | is.factor(FBrOrd[[keepVar]])){
      Bxplt(
        x = FBrOrd[[keepVar]],
        y =  FBrOrd$FBratio,
        fill = FBrOrd[match(unique(FBrOrd[[keepVar]]), FBrOrd[[keepVar]]), paste(keepVar, "col", sep = "_")],
        colpts = FBrOrd[[paste(keepVar, "col", sep = "_")]],
        las = ifelse(length(unique(FBrOrd[[keepVar]])) > 4, 2, 1),
        boxwex = 0.5,
        ylab = "Firmicutes/Bacteroidota ratio",
        cex.axis = ifelse(length(unique(FBrOrd[[keepVar]])) > 4, 0.7, 1.1),
        main = paste(i, keepVar, sep = " : "),
        yscale = 6)
    }else if(is.numeric(FBrOrd[[j]])) {
      plot(x = FBrOrd[[keepVar]],
           y =  FBrOrd$FBratio,
           xlab = keepVar,
           ylab = "Firmicutes/Bacteroidota ratio",
           pch = 19,
           main = paste(i, keepVar, sep = " : "))
      abline(coef = c( 
        BestMod$ModSum$coefficients["(Intercept)", "Estimate"],
        BestMod$ModSum$coefficients[grep(keepVar, rownames(BestMod$ModSum$coefficients)), "Estimate"]), 
        col = "firebrick")
    }
  }
  dev.off()

}

# save statistical test outputs
saveRDS(RatioRes, file = file.path(DirSave, "Signif_Effects", "LMMResTaxRatio_FB.rds"), compress = TRUE)

# test for BP ratio
RatioRes <- list()

for (i in names(normalizedDat)) {
  
  # Convert to phyloseq object 
  ## makes otu table for phyloseq
  read.phylo <-
    phyloseq::otu_table(normalizedDat[[i]], taxa_are_rows = T)
  
  ## make sample table for phyloseq
  sample.phylo <- phyloseq::sample_data(labFishPhy[["samples"]])
  
  ## merge phyloseq objects
  PhyloDat  <- 
    phyloseq::phyloseq(read.phylo, sample.phylo)
  
  # extract the count for selected phylums
  Bacteroidota <- PhyloDat@otu_table[which(rownames(PhyloDat@otu_table) == "Bacteroidota")]
  Proteobacteria <- PhyloDat@otu_table[which(rownames(PhyloDat@otu_table) == "Proteobacteria")]
  
  ##############################################
  # compute Bacteroidota/Proteobacteria ratio  #
  ##############################################
  PBr <- data.frame(Sample = colnames(Proteobacteria), Proteobacteria = as.vector(Proteobacteria), Bacteroidota = as.vector(Bacteroidota))
  PBr[["PBratio"]] <- PBr[["Bacteroidota"]]/PBr[["Proteobacteria"]]
  PBr <- cbind(PBr, as.data.frame(PhyloDat@sam_data))
  PBr$PBratio[which(is.infinite(PBr$PBratio))] <- NA
  PBr <- PBr[-which(is.na(PBr$PBratio)),]
  
  # test the effect of treatment on the ratio
  PBr[["Bac"]] <- as.character(PBr[["Bac"]])
  BestMod <- DredgedLMM(PBr, 
                        expVar = "PBratio",
                        respExpr = "(Contam + ContPop + Inj)^3 + Taille_mm_M",
                        random = "(1 |Session / Bac) + (1 |Pop)",
                        method = "weight",
                        rank = "AICc")
  
  
  DirSave <- file.path(savingDir, "Normalized_Data")
  if (length(list.dirs(file.path(DirSave, "Signif_Effects"))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects"))
  }
  if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio"))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio"))
  }
  if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio", i))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio", i))
  }
  if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio", i, "PBratio"))) == 0) {
    dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio", i, "PBratio"))
  }
  toSave <- file.path(DirSave, "Signif_Effects", "TaxRatio", i, "PBratio", paste0("Mod_Perf_", "PBratio","_",  i, ".tif" ))
  if(file.exists(toSave)){
    unlink(toSave)
  } 
  
  grDevices::tiff(filename = toSave,
                  width = 1200,
                  height = 1200,
                  res = 100,
                  compression = "lzw",
                  pointsize = 8)
  print(performance::check_model(BestMod[["ModL"]]))
  dev.off()
  
  # store the result of the model
  RatioRes[[i]][["PBratio"]] <- BestMod
  
  # display the significant results 
  toKeep <- rownames(BestMod$ModAnov)[-1]
  
  if(length(toKeep) == 0){
    next
  }
  ##### save the figs representing difference of PCOA Scores (composition among significant treatments)
  toSave <- file.path(DirSave, "Signif_Effects", "TaxRatio", i, "PBratio", paste0("SignifPlot", "PBratio","_",  i, ".tif" ))
  if(file.exists(toSave)){
    unlink(toSave)
  } 
  grDevices::tiff(filename = toSave,
                  width = 1200,
                  height = 1200,
                  res = 300,
                  compression = "lzw",
                  pointsize = 8)
  par(mfrow = c(ifelse(ceiling(sqrt(length(toKeep)))^2 > length(toKeep), ceiling(sqrt(length(toKeep))),2), ceiling(sqrt(length(toKeep)))), 
      mar = c(4, 4, 4, 2))
  for(j in seq_along(toKeep)){ 
    if(!length(grep(":", toKeep[[j]])) > 0){
      keepVar <- toKeep[[j]]
    }else{
      keepVar <-  gsub(":", ".", toKeep[[j]])
    }
    PBrOrd <- PBr[order(PBr[[keepVar]]), ]
    if(is.character(PBrOrd[[keepVar]]) | is.factor(PBrOrd[[keepVar]])){
      Bxplt(
        x = PBrOrd[[keepVar]],
        y =  PBrOrd$PBratio,
        fill = PBrOrd[match(unique(PBrOrd[[keepVar]]), PBrOrd[[keepVar]]), paste(keepVar, "col", sep = "_")],
        colpts = PBrOrd[[paste(keepVar, "col", sep = "_")]],
        las = ifelse(length(unique(PBrOrd[[keepVar]])) > 4, 2, 1),
        boxwex = 0.5,
        ylab = "Proteobacteria/Bacteroidota ratio",
        cex.axis = ifelse(length(unique(PBrOrd[[keepVar]])) > 4, 0.7, 1.1),
        main = paste(i, keepVar, sep = " : "),
        yscale = 6)
    }else if(is.numeric(PBrOrd[[j]])) {
      plot(x = PBrOrd[[keepVar]],
           y =  PBrOrd$PBratio,
           xlab = keepVar,
           ylab = "Proteobacteria/Bacteroidota ratio",
           pch = 19,
           main = paste(i, keepVar, sep = " : "))
      abline(coef = c( 
        BestMod$ModSum$coefficients["(Intercept)", "Estimate"],
        BestMod$ModSum$coefficients[grep(keepVar, rownames(BestMod$ModSum$coefficients)), "Estimate"]), 
        col = "firebrick")
    }
  }
  dev.off()
  
}

# save statistical test outputs
saveRDS(RatioRes, file = file.path(DirSave, "Signif_Effects", "LMMResTaxRatio_BP.rds"), compress = TRUE)


###############################
#                             #                     
#      Summary table          #
#                             #
###############################

# for FB ratio
res <- readRDS(file = file.path(DirSave, "Signif_Effects", "LMMResTaxRatio_FB.rds"))

# create a temporary table grouping the names of the best models and the name of the corresponding response variable
BestMods <- data.frame(
  mods = names(res),
  modsId = c(
    "Original",
    "Original (log-transformed)",
    "Proportion",
    "Proportion (log-transformed)",
    "CSS",
    "CSS (log-transformed)",
    #"Rarefied",
    #"Rarefied (log-transformed)",
    "DESEQ",
    "DESEQ (log-transformed)",
    "TMM",
    "TMM (log-transformed)"
  )
)

# retrieve the sample size and the Rsquared
i = "FBRatio"
  BestMods[[paste("n", i, sep = "_")]] <- rep(NA, nrow(BestMods))
  BestMods[[paste("R2m", i, sep = "_")]] <- rep(NA, nrow(BestMods))
  BestMods[[paste("R2c", i, sep = "_")]] <- rep(NA, nrow(BestMods))
  for (r in seq(nrow(BestMods))) {
    # n
    BestMods[r, paste("n", i, sep = "_")] <-
      paste0("n = ", res[[BestMods[r, "mods"]]][[i]]$ModSum$devcomp$dims[["n"]])
    # Rsquared
    BestMods[r, paste("R2m", i, sep = "_")] <-
      paste0("R2m = ", signif(res[[BestMods[r, "mods"]]][[i]]$Rsquared[[1]], digits = 3))
    BestMods[r, paste("R2c", i, sep = "_")] <-
      paste0("R2c = ", signif(res[[BestMods[r, "mods"]]][[i]]$Rsquared[[2]], digits = 3))
  }

# retrieve the useful informations from the models
## create a quick function to retrieve the results for factors with more than two levels 
find_matches <- function(strings) {
  # Extract the first three characters of each string
  prefixes <- substr(strings, 1, 3)
  # Find unique prefixes (as each unique prefix represents a different starting sequence)
  unique_prefixes <- unique(prefixes)
  # Check the frequency of each unique prefix
  prefix_count <- table(prefixes)
  # Find out if any of the unique prefixes occurs more than once (indicating a common sequence)
  common_prefixes <- unique_prefixes[prefix_count[unique_prefixes] > 1]
  # If there are common prefixes, return the strings that have them; otherwise, return a message
  if (length(common_prefixes) > 0) {
    matching_strings <- lapply(common_prefixes, function(prefix) {
      strings[grep(paste0("^", prefix), strings)]
    })
    return(matching_strings)
  } else {
    return("NA")
  }
}

# retrieve information from model summary and anova table
options(scipen = 999)
  ValuesRes <- data.frame()
  for (r in seq(nrow(BestMods))) {
    coef <- res[[BestMods[r, "mods"]]][[i]]$ModSum[["coefficients"]]
    df <- res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Df"]]
    Chisq = res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Chisq"]]
    p.value = res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Pr(>Chisq)"]]
    if(nrow(coef) != length(df)){
      match <- unlist(find_matches(rownames(coef)))
      # deal with interactions
      isCont <- grep("Cont", match)
      if(length(isCont)>0 & length(grep(":", isCont)) == 0 ){
        match <- match[-isCont]
      }
      char <- substr(match[1], 1, 3)
      index <- min(grep(char, rownames(coef)))
      df <- append(df, df[index], after=index)
      Chisq <- append(Chisq, Chisq[index], after=index)
      p.value <- append(p.value, p.value[index], after=index)
    }
    Values <- cbind(coef, df, Chisq, p.value)
    Values <- signif(Values, digits = 3)
    Values <- cbind(rownames(Values), Values)
    rownames(Values) <- NULL
    Values <- cbind(rep(BestMods[r, "mods"], nrow(Values)), Values)
    ValuesRes <- rbind(ValuesRes, Values)
  }
  # leave the two first column without names
  names(ValuesRes)[c(1, 2)] <- c("", "")
  
  # specify the new predictor names (to display in the table)
  PredictorsNames <- data.frame(toreplace = unique(ValuesRes[[2]]), replacement = NA)
  PredictorsNames$replacement[PredictorsNames$toreplace == "(Intercept)"] <- "Intercept"
  PredictorsNames$replacement[PredictorsNames$toreplace == "InjPBS"] <- "Imm. Chall. (PBS)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "SexeI"] <- "Sex (I)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "SexeM"] <- "Sex (M)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContamNC"] <- "Contamination (NC)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContPopNC"] <- "Origin (NC)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContamNC:ContPopNC"] <- "Contam. (NC):Origin (NC)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContamNC:InjPBS"] <- "Contam. (NC):Imm. Chall. (PBS)"
  PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
  
  # rename the predictors
  for (pn in PredictorsNames[["toreplace"]]) {
    ValuesRes[[2]][which(ValuesRes[[2]] == pn)] <-
      PredictorsNames[["replacement"]][which(PredictorsNames[["toreplace"]] == pn)]
  }
  
  # replace p.value with large decimal by approximations
  ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] < 0.0001)] <-
    "<0.0001"
  ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] > 0.0001 &
                                 ValuesRes[["p.value"]] < 0.01)] <- "<0.01"
  ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] > 0.01 &
                                 ValuesRes[["p.value"]] < 0.05)] <- "<0.05"

options(scipen = 0)
str(ValuesRes)

# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalization methods)
## create a directory to store the table containing LMM results for alpha diversity indices according to selected normalization
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio", "LMMFBRatioGutBestMod"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio", "LMMFBRatioGutBestMod"))
} 

  l <- vector("list", nrow(BestMods))
  l[[1]] <- kableExtra::kable_classic_2(
    kableExtra::kbl(ValuesRes[-1], 
                    align = "c", 
                    caption = paste0("<span style='font-size:12px; font-weight: bold; font-style: italic'>", "FBRatio", "</span>")),
    bootstrap_options = c("striped", "hover", "condensed", "responsive"),
    full_width = F,
    html_font = "arial",
    font_size = 10
  )

  for (i in seq(nrow(BestMods))) {
    l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                        if (is.na(BestMods[i, paste("R2m", "FBRatio", sep = "_")]) &
                                            is.na(BestMods[i, paste("R2c", "FBRatio", sep = "_")]) &
                                            is.na(BestMods[i, paste("n", "FBRatio", sep = "_")])) {
                                          paste(BestMods[i, "modsId"], BestMods[i, paste("n", "FBRatio", sep = "_")], sep = " | ")
                                        } else{
                                          paste(BestMods[i, "modsId"], BestMods[i, paste("n", "FBRatio", sep = "_")], BestMods[i, paste("R2m", "FBRatio", sep = "_")], BestMods[i, paste("R2c", "FBRatio", sep = "_")], sep = " | ")
                                        },
                                        min(which(ValuesRes ==  BestMods[i, "mods"])),
                                        max(which(ValuesRes ==  BestMods[i, "mods"])),
                                        label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
  }
  l[[length(l)]]
  kableExtra::save_kable(l[[length(l)]], file.path(DirSave, "Signif_Effects", "TaxRatio", "LMMFBRatioGutBestMod", paste0(paste("LMMFBRatioGutBestMod", "FBRatio", sep = "_"), ".HTML")))



# for BP ratio
res <- readRDS(file = file.path(DirSave, "Signif_Effects", "LMMResTaxRatio_BP.rds"))

# create a temporary table grouping the names of the best models and the name of the corresponding response variable
BestMods <- data.frame(
  mods = names(res),
  modsId = c(
    "Original",
    "Original (log-transformed)",
    "Proportion",
    "Proportion (log-transformed)",
    "CSS",
    "CSS (log-transformed)",
    #"Rarefied",
    #"Rarefied (log-transformed)",
    "DESEQ",
    "DESEQ (log-transformed)",
    "TMM",
    "TMM (log-transformed)"
  )
)

# retrieve the sample size and the Rsquared
i = "PBratio"
  BestMods[[paste("n", i, sep = "_")]] <- rep(NA, nrow(BestMods))
  BestMods[[paste("R2m", i, sep = "_")]] <- rep(NA, nrow(BestMods))
  BestMods[[paste("R2c", i, sep = "_")]] <- rep(NA, nrow(BestMods))
  for (r in seq(nrow(BestMods))) {
    # n
    BestMods[r, paste("n", i, sep = "_")] <-
      paste0("n = ", res[[BestMods[r, "mods"]]][[i]]$ModSum$devcomp$dims[["n"]])
    # Rsquared
    BestMods[r, paste("R2m", i, sep = "_")] <-
      paste0("R2m = ", signif(res[[BestMods[r, "mods"]]][[i]]$Rsquared[[1]], digits = 3))
    BestMods[r, paste("R2c", i, sep = "_")] <-
      paste0("R2c = ", signif(res[[BestMods[r, "mods"]]][[i]]$Rsquared[[2]], digits = 3))
  }

# retrieve the useful informations from the models
## create a quick function to retrieve the results for factors with more than two levels 
find_matches <- function(strings) {
  # Extract the first three characters of each string
  prefixes <- substr(strings, 1, 3)
  # Find unique prefixes (as each unique prefix represents a different starting sequence)
  unique_prefixes <- unique(prefixes)
  # Check the frequency of each unique prefix
  prefix_count <- table(prefixes)
  # Find out if any of the unique prefixes occurs more than once (indicating a common sequence)
  common_prefixes <- unique_prefixes[prefix_count[unique_prefixes] > 1]
  # If there are common prefixes, return the strings that have them; otherwise, return a message
  if (length(common_prefixes) > 0) {
    matching_strings <- lapply(common_prefixes, function(prefix) {
      strings[grep(paste0("^", prefix), strings)]
    })
    return(matching_strings)
  } else {
    return("NA")
  }
}

# retrieve information from model summary and anova table
options(scipen = 999)
  ValuesRes <- data.frame()
  for (r in seq(nrow(BestMods))) {
    coef <- res[[BestMods[r, "mods"]]][[i]]$ModSum[["coefficients"]]
    df <- res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Df"]]
    Chisq = res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Chisq"]]
    p.value = res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Pr(>Chisq)"]]
    if(nrow(coef) != length(df)){
      match <- unlist(find_matches(rownames(coef)))
      # deal with interactions
      isCont <- grep("Cont", match)
      if(length(isCont)>0 & length(grep(":", isCont)) == 0 ){
        match <- match[-isCont]
      }
      char <- substr(match[1], 1, 3)
      index <- min(grep(char, rownames(coef)))
      df <- append(df, df[index], after=index)
      Chisq <- append(Chisq, Chisq[index], after=index)
      p.value <- append(p.value, p.value[index], after=index)
    }
    Values <- cbind(coef, df, Chisq, p.value)
    Values <- signif(Values, digits = 3)
    Values <- cbind(rownames(Values), Values)
    rownames(Values) <- NULL
    Values <- cbind(rep(BestMods[r, "mods"], nrow(Values)), Values)
    ValuesRes <- rbind(ValuesRes, Values)
  }
  # leave the two first column without names
  names(ValuesRes)[c(1, 2)] <- c("", "")
  
  # specify the new predictor names (to display in the table)
  PredictorsNames <- data.frame(toreplace = unique(ValuesRes[[2]]), replacement = NA)
  PredictorsNames$replacement[PredictorsNames$toreplace == "(Intercept)"] <- "Intercept"
  PredictorsNames$replacement[PredictorsNames$toreplace == "InjPBS"] <- "Imm. Chall. (PBS)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "SexeI"] <- "Sex (I)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "SexeM"] <- "Sex (M)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContamNC"] <- "Contamination (NC)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContPopNC"] <- "Origin (NC)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContamNC:ContPopNC"] <- "Contam. (NC):Origin (NC)"
  PredictorsNames$replacement[PredictorsNames$toreplace == "ContamNC:InjPBS"] <- "Contam. (NC):Imm. Chall. (PBS)"
  PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
  
  # rename the predictors
  for (pn in PredictorsNames[["toreplace"]]) {
    ValuesRes[[2]][which(ValuesRes[[2]] == pn)] <-
      PredictorsNames[["replacement"]][which(PredictorsNames[["toreplace"]] == pn)]
  }
  
  # replace p.value with large decimal by approximations
  ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] < 0.0001)] <-
    "<0.0001"
  ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] > 0.0001 &
                                 ValuesRes[["p.value"]] < 0.01)] <- "<0.01"
  ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] > 0.01 &
                                 ValuesRes[["p.value"]] < 0.05)] <- "<0.05"

options(scipen = 0)
str(ValuesRes)

# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalization methods)
## create a directory to store the table containing LMM results for alpha diversity indices according to selected normalization
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "TaxRatio", "LMMPBratioGutBestMod"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "TaxRatio", "LMMPBratioGutBestMod"))
} 

l <- vector("list", nrow(BestMods))
l[[1]] <- kableExtra::kable_classic_2(
  kableExtra::kbl(ValuesRes[-1], 
                  align = "c", 
                  caption = paste0("<span style='font-size:12px; font-weight: bold; font-style: italic'>", "PBratio", "</span>")),
  bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  full_width = F,
  html_font = "arial",
  font_size = 10
)

for (i in seq(nrow(BestMods))) {
  l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                      if (is.na(BestMods[i, paste("R2m", "PBratio", sep = "_")]) &
                                          is.na(BestMods[i, paste("R2c", "PBratio", sep = "_")]) &
                                          is.na(BestMods[i, paste("n", "PBratio", sep = "_")])) {
                                        paste(BestMods[i, "modsId"], BestMods[i, paste("n", "PBratio", sep = "_")], sep = " | ")
                                      } else{
                                        paste(BestMods[i, "modsId"], BestMods[i, paste("n", "PBratio", sep = "_")], BestMods[i, paste("R2m", "PBratio", sep = "_")], BestMods[i, paste("R2c", "PBratio", sep = "_")], sep = " | ")
                                      },
                                      min(which(ValuesRes ==  BestMods[i, "mods"])),
                                      max(which(ValuesRes ==  BestMods[i, "mods"])),
                                      label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
}
l[[length(l)]]
kableExtra::save_kable(l[[length(l)]], file.path(DirSave, "Signif_Effects", "TaxRatio", "LMMPBRatioGutBestMod", paste0(paste("LMMPBRatioGutBestMod", "PBratio", sep = "_"), ".HTML")))
