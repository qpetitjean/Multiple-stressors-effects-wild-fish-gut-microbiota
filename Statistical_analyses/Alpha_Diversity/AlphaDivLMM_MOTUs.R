# Script Title: Compute Alpha Diversity Indices and Model Treatment Effects in Fish Gut Microbiota
#      
# Author: Quentin PETITJEAN
# Date Created: 03/2023
# Last Modified: 05/05/2024
# ==============================================================================
# Requirements: 
# - R version 4.2.3
# - Packages:
#   - metabaR v1.0.0: For processing and subsetting metabarcoding datasets.
#   - ape v5.7-1:: For phylogenetic tree manipulation.
#   - car v3.1-2: For computing type-III ANOVA on linear models.
#   - kableExtra v1.3.4: For creating enhanced HTML tables of statistical results.
# - Source Files:
#   - multiNorm.R: Function for applying multiple normalization methods.
#   - DredgedLMM.R: Function to perform model dredging for LMM and retrieve summaries.
#   - alphaDiv.R: Function to compute alpha diversity indices.
#   - BxpltFunc.R: Function to generate customized boxplots.
# ==============================================================================
# Script Overview:
# This script processes fish gut microbiota data by integrating experimental metadata,
# a phylogenetic tree, and a cleaned metabar object. The workflow includes:
# 1. Importing custom functions via GitHub.
# 2. Loading the full experimental dataset, a phylogenetic tree, and the fish gut metabar dataset.
# 3. Computing multiple normalization methods on the fish gut data.
# 4. Calculating alpha diversity indices (e.g., Richness, Shannon, Simpson, PD) based on the
#    normalized data and phylogenetic tree when needed.
# 5. Running linear mixed models (LMM) using a dredged model selection approach to test treatment effects,
#    incorporating factors such as Contamination, Population, Injection, and Sex.
# 6. Generating diagnostic plots and saving model performance outputs.
# 7. Compiling HTML summary tables for LMM results across normalization methods.
# ==============================================================================
# Usage:
# 1. Update the 'savingDir' variable with the correct path for your input/output files.
# 2. Ensure that all required data files (CSV, RDS, NWK) and source URLs for custom functions are accessible.
# 3. Install the necessary R packages if not already installed.
# 4. Run the script in an R environment to perform normalization, compute alpha diversity indices,
#    run LMM analyses, and generate summary outputs and diagnostic plots.
# 5. The script outputs include:
#    - An RDS file containing LMM results.
#    - Diagnostic plot files (TIF format) for model performance.
#    - HTML summary tables for each normalization method.
#    - Boxplots visualizing treatment effects on alpha diversity indices.
# ==============================================================================

##############################################
#       	Install needed packages            #
##############################################

if(!require(metabaR)){
  install.packages("metabaR")
}
if(!require(ape)){
  install.packages("ape")
}
if(!require(car)){
  install.packages("car")
}
if(!require(kableExtra)){
  install.packages("kableExtra")
}

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported and saved
savingDir <- "w:/POSTDOC_INP_GOLFECH_2023/Outputs"

##############################################
#       	Import some custom functions       #
##############################################
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/multiNorm.R") # import a function used to compute multiple normalization methods
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/DredgedLMM.R")  # import a function used to retrieve the summary, anova table and Rsquared of the most straightforward mixed effect model using rank selection (e.g., AICc)
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/alphaDiv.R") # import a function compute alpha diversity indices from a matrix containing read counts per sample
source("https://raw.githubusercontent.com/qpetitjean/Multiple-stressors-effects-wild-fish-gut-microbiota/main/R_Func/BxpltFunc.R") # import a function compute draw customized boxplots

#######################
# Import the data     #
#######################

# import the full dataset from the experiment (treatment, sex and other fish variable)
FullDat <-
  read.csv2(file.path(savingDir, "Data/DesignData", "DesignData.csv"), dec = ".", sep = ";")

# import the phylogenetic tree (.nwk)
PhyloTree <-
  ape::read.tree(file.path(savingDir, "Data/PhyloTree/PhyloTree.nwk"))

# import the cleaned dataset (metabaR object) and select only gut samples
labFish <- readRDS(file.path(savingDir, "Data/CleanedData", "fguts_Bact_agg_MergedRep.RDS"))
labFish <- metabaR::subset_metabarlist(labFish,
                                       table = "pcrs",
                                       indices = labFish$pcrs$matrix == "fish_gut")
rownames(labFish[["reads"]]) <- labFish$samples$Num_prlvt_Euth

#######################################
# compute multiple normalization      #
#######################################
normalizedDat <- multiNorm(
  dataList = labFish,
  norm = c(
    "original",
    "original.log",
    "read.prop",
    "read.prop.log",
    "read.CSS",
    "read.CSS.log",
    "read.rare",
    "read.rare.log",
    "read.deseq",
    "read.deseq.log",
    "read.tmm",
    "read.tmm.log"
  )
)

#############################################################
# compute alpha div. indices over the normalized data       #
#############################################################
indices = c(
  "Richness",
  "Shannon",
  "Shannon.exp",
  "Simpson",
  "PD"
)

alphaDivindices <- lapply(normalizedDat, function(x)
  alphaDiv(
    Data = x,
    indices = indices,
    PhyloTree = PhyloTree
  )
)

################################################################
# Run LLM on the normalized dataset and over the indices       #
################################################################
# create a directory and subdirectories to save the output plots and data
if (length(list.dirs(file.path(savingDir, "Normalized_Data"))) == 0) {
  dir.create(file.path(savingDir, "Normalized_Data"))
}
DirSave <- file.path(savingDir, "Normalized_Data")
if (length(list.dirs(file.path(DirSave, "Signif_Effects"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects"))
}
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "LMM"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "LMM"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "LMM", "AlphaDivPlots"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "LMM", "AlphaDivPlots"))
} 

# create several versions of the full dataset by merging the alphadiv indices to the full dataset of the experiment
FinalDat <- lapply(alphaDivindices, function(x)
  merge(
    x = FullDat,
    y = x,
    by.x = "Ind",
    by.y = "Samples",
    all = TRUE
  ))

res <- list()
# iterate trough the various dataset (normalization) to test the effect of the experimental design on alpha diversity indices
for (i in seq_along(FinalDat)){
  Temp <- FinalDat[[i]]
  for (j in seq_along(indices)){
    if(length(which(is.na(Temp[indices[[j]]]))) > 0 && 
       nrow(Temp[-which(is.na(Temp[indices[[j]]])), ]) == 0){
      next
    }
    # dredge the full model and retrieve the summary, anova table and R2 of the most straightforward one 
    res[[names(FinalDat)[i]]][[indices[j]]] <- DredgedLMM(FinalDat[[i]], 
                                                          expVar = indices[[j]],
                                                          respExpr = "(Contam + ContPop + Inj)^3 + Taille_mm_M + Sexe",
                                                          random = "(1 |Session / Bac) + (1 | Pop)",
                                                          method = "weight",
                                                          rank = "AICc")
    # save performance plot (LMM diagnostic)
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "LMM", "AlphaDivPlots", names(FinalDat)[i]))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "LMM", "AlphaDivPlots", names(FinalDat)[i]))
    } 
    toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots", "LMM", "AlphaDivPlots", names(FinalDat)[i], paste0(paste("Mod_Perf_AlphaDivLMM", names(FinalDat)[i], indices[j], sep = "_"), ".tif" ))
    if(file.exists(toSave)){
      unlink(toSave)
    }
    grDevices::tiff(filename = toSave,
                    width = 1200,
                    height = 1200,
                    res = 100,
                    compression = "lzw",
                    pointsize = 8)
    print(performance::check_model(res[[names(FinalDat)[i]]][[indices[j]]]$ModL))
    dev.off()
  }
}

# save statistical test outputs
saveRDS(res, file = file.path(DirSave, "Signif_Effects", "LMMAlphaDivRes.rds"))

## plot the results
for (i in seq_along(FinalDat)){
  Temp <- FinalDat[[i]]
  for(j in seq_along(indices)){
    if(length(which(is.na(Temp[indices[[j]]]))) > 0 && 
       nrow(Temp[-which(is.na(Temp[indices[[j]]])), ]) == 0){
      next
    }
    toSave <- file.path(DirSave, "Signif_Effects/Signif_Plots/LMM/AlphaDivPlots", names(FinalDat)[i], gsub(":", "X", paste0(paste("AlphaDiv", names(FinalDat)[i], gsub("\\.", "_", indices[[j]]), sep = "_"), ".tif")))
    if(file.exists(toSave)){
    unlink(toSave)
    }
    grDevices::tiff(filename = toSave,
                width = 1200,
                height = 1200,
                res = 300,
                compression = "lzw",
                pointsize = 8)
par(mfrow=c(2,2))

##### plot differences among sex
Bxplt(
  x = Temp[["Sexe"]],
  y = Temp[[indices[[j]]]],
  fill = c(
    adjustcolor("#7B68EE", alpha.f = 0.3),
    adjustcolor("#B0C4DE", alpha.f = 0.3),
    adjustcolor("#87CEEB", alpha.f = 0.3)
  ),
  colpts = Temp[["Sexe_col"]],
  boxwex = 0.5,
  ylab = indices[[j]],
  cex.axis = 0.9,
  ymin = 0,
  yscale = 6
)

##### plot differences among Contamination treatment
newOrder <- c("NC", "C")
Temp[["ContamOrd"]] <- factor(Temp[["Contam"]], levels = newOrder)

Bxplt(
  x = Temp[["ContamOrd"]],
  y = Temp[[indices[[j]]]],
  fill = c(
    adjustcolor("#006400", alpha.f = 0.3),
    adjustcolor("#B22222", alpha.f = 0.3)
  ),
  colpts = Temp[["Contam_col"]],
  boxwex = 0.5,
  ylab = indices[[j]],
  cex.axis = 0.9,
  ymin = 0,
  yscale = 6
)

##### plot differences among Injection treatment
newOrder <- c("PBS", "AMIX")
Temp[["InjOrd"]] <- factor(Temp[["Inj"]], levels = newOrder)

Bxplt(
  x = Temp[["InjOrd"]],
  y = Temp[[indices[[j]]]],
  fill = c(
    adjustcolor("#808080", alpha.f = 0.3),
    adjustcolor("#470C8E", alpha.f = 0.3)
  ),
  colpts = Temp[["Inj_col"]],
  boxwex = 0.5,
  ylab = indices[[j]],
  cex.axis = 0.9,
  ymin = 0,
  yscale = 6
)

##### plot differences among populations
newOrder <- c("ELEV", "ARIMAS", "CELCAB", "AUSCOR", "RIOU")
Temp[["PopOrd"]] <- factor(Temp[["Pop"]], levels = newOrder)

Bxplt(
  x = Temp[["PopOrd"]],
  y = Temp[[indices[[j]]]],
  fill = unique(Temp[["Pop_col"]])[match(newOrder, unique(Temp[["Pop"]]))],
  colpts = Temp[["Pop_col"]],
  boxwex = 0.5,
  ylab = indices[[j]],
  cex.axis = 0.9,
  las = 2,
  ymin = 0,
  yscale = 6
)
dev.off()
}
}

res <- readRDS(file.path(DirSave, "Signif_Effects", "LMMAlphaDivRes.rds"))

################################################################
# summarize the results in html table for each normalization  #
################################################################

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
    "Rarefied",
    "Rarefied (log-transformed)",
    "DESEQ",
    "DESEQ (log-transformed)",
    "TMM",
    "TMM (log-transformed)"
  )
)

# retrieve the sample size and the Rsquared
for(i in indices){
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
Val <- setNames(lapply(indices, function(i) {
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
return(ValuesRes)
}), indices)
options(scipen = 0)

# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalization methods)
## create a directory to store the table containing LMM results for alpha diversity indices according to selected normalization
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "LMMAlphaDivBestModTab"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "LMMAlphaDivBestModTab"))
} 
for(m in indices){
  ValIndex <- which(names(Val) == m)
l <- vector("list", nrow(BestMods))
l[[1]] <- kableExtra::kable_classic_2(
  kableExtra::kbl(Val[[ValIndex]][-1], 
                  align = "c", 
                  caption = paste0("<span style='font-size:12px; font-weight: bold; font-style: italic'>", m, "</span>")),
  bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  full_width = F,
  html_font = "arial",
  font_size = 10
)
for (i in seq(nrow(BestMods))) {
  l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                      if (is.na(BestMods[i, paste("R2m", m, sep = "_")]) &
                                          is.na(BestMods[i, paste("R2c", m, sep = "_")])) {
                                        paste(BestMods[i, "modsId"], BestMods[i, paste("n", m, sep = "_")], sep = " | ")
                                      } else{
                                        paste(BestMods[i, "modsId"], BestMods[i, paste("n", m, sep = "_")], BestMods[i, paste("R2m", m, sep = "_")], BestMods[i, paste("R2c", m, sep = "_")], sep = " | ")
                                      },
                                      min(which(Val[[ValIndex]][[1]] ==  BestMods[i, "mods"])),
                                      max(which(Val[[ValIndex]][[1]] ==  BestMods[i, "mods"])),
                                      label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
}
l[[length(l)]]
kableExtra::save_kable(l[[length(l)]], file.path(DirSave, "Signif_Effects", "LMMAlphaDivBestModTab", paste0(paste("LMMAlphaDivBestMod", m, sep = "_"), ".HTML")))
}

# generate the table including all metric for the original data log transformed (Table 2)
## retrieve the summary from the indices of the original Log transformed data
dataNorm <- "original.log"

BestModPart <- BestMods[which(BestMods$mods == dataNorm),]
BestModPart <- BestModPart[,-c(1,2)]
BestModPart <- setNames(split.default(BestModPart, gl(ncol(BestModPart) / 3, 3)), indices)

OrLog <- lapply(Val, function(x)
  x[which(x[[1]] == dataNorm), ])
OrLogDf <- do.call(rbind,OrLog)
OrLogDf$index <- gsub("\\.\\d+", "", rownames(OrLogDf))
rownames(OrLogDf) <- NULL

l <- vector("list", length(BestModPart))
l[[1]] <- kableExtra::kable_classic_2(
  kableExtra::kbl(
    OrLogDf[-c(1, length(OrLogDf))],
    align = "c",
    caption = paste0(
      "<span style='font-size:12px; font-weight: bold; font-style: italic'>",
      dataNorm,
      "</span>"
    )
  ),
  bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  full_width = F,
  html_font = "arial",
  font_size = 10
)

for (i in seq_along(BestModPart)) {
  l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                      if (is.na(BestModPart[[i]][, paste("R2m", names(BestModPart)[[i]], sep = "_")]) &
                                          is.na(BestModPart[[i]][, paste("R2c", names(BestModPart)[[i]], sep = "_")])) {
                                        paste(names(BestModPart)[[i]], BestModPart[[i]][, paste("n", names(BestModPart)[[i]], sep = "_")], sep = " | ")
                                      } else{
                                        paste(names(BestModPart)[[i]], BestModPart[[i]][, paste("n", names(BestModPart)[[i]], sep = "_")], BestModPart[[i]][, paste("R2m", names(BestModPart)[[i]], sep = "_")], BestModPart[[i]][, paste("R2c", names(BestModPart)[[i]], sep = "_")], sep = " | ")
                                      },
                                      min(which(OrLogDf[["index"]] ==  names(BestModPart)[[i]])),
                                      max(which(OrLogDf[["index"]] ==  names(BestModPart)[[i]])),
                                      label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
}

l[[length(l)]]
kableExtra::save_kable(l[[length(l)]], file.path(DirSave, "Signif_Effects", "LMMAlphaDivBestModTab", paste0(
  paste("LMMAlphaDivBestMod", dataNorm, sep = "_"), ".HTML"
)))


# compute mean between contamination treatment for all alpha diversity indices (reported in the results section)
MeanAlphaDivContam <- setNames(lapply(indices, function(y)
  aggregate(
    FinalDat$original.log[[y]],
    list(FinalDat$original.log$Contam),
    FUN = function(x) {
      mean(x, na.rm = T)
    }
  )), indices)

MeanAlphaDivContam

# compute se between contamination treatment for all alpha diversity indices (reported in the results section)
SeAlphaDivContam <- setNames(lapply(indices, function(y)
    stats::aggregate(
      FinalDat$original.log[[y]],
      list(FinalDat$original.log$Contam),
      FUN = function(x) {
        sd(x, na.rm = T)
      }
    )$x / stats::aggregate(FinalDat$original.log[[y]], list(FinalDat$original.log$Contam), function(x)
      sqrt(sum(!is.na(x))))$x), indices)

SeAlphaDivContam

# display the boxplot of richeness differences among Contamination treatment for the original log transformed data (Figure 1A)
# Reorder the factor levels
newOrder <- c("NC", "C")
FinalDat$original.log[["ContamOrd"]] <- factor(FinalDat$original.log[["Contam"]], levels = newOrder)
toSave <- file.path(DirSave, "Signif_Effects/Signif_Plots/LMM/AlphaDivPlots/BP_RICH_LogORDat.tif")
if(file.exists(toSave)){
  unlink(toSave)
}
grDevices::tiff(filename = toSave,
                width = 1200,
                height = 1200,
                res = 300,
                compression = "lzw",
                pointsize = 11)
Bxplt(
  x = FinalDat$original.log[["ContamOrd"]],
  y = FinalDat$original.log[["Richness"]],
  fill = c(
    adjustcolor("#006400", alpha.f = 0.3),
    adjustcolor("#B22222", alpha.f = 0.3)
  ),
  colpts = FinalDat$original.log[["Contam_col"]],
  boxwex = 0.5,
  ylab = "Richness (# of MOTUs)",
  cex.axis = 0.9,
  ymin = 0,
  yscale = 6
)
dev.off()


# display the boxplot of differences among Contamination treatment for the original log transformed data for all indices
# Reorder the factor levels
newOrder <- c("NC", "C")
FinalDat$original.log[["ContamOrd"]] <- factor(FinalDat$original.log[["Contam"]], levels = newOrder)
toSave <- file.path(DirSave, "Signif_Effects/Signif_Plots/LMM/AlphaDivPlots/BP_AllIndices_LogORDat.tif")
if(file.exists(toSave)){
  unlink(toSave)
}
grDevices::tiff(filename = toSave,
                width = 1200*5,
                height = 1200,
                res = 300,
                compression = "lzw",
                pointsize = 22)
par(mfrow=c(1,5))
for(i in seq_along(indices)){
ifelse(rep(i == 1, 2), par(mar=c(5, 4, 4, 2) + 0.1), par(mar=c(4, 2, 4, 2)))  
Bxplt(
  x = FinalDat$original.log[["ContamOrd"]],
  y = FinalDat$original.log[[indices[i]]],
  fill = c(
    adjustcolor("#006400", alpha.f = 0.3),
    adjustcolor("#B22222", alpha.f = 0.3)
  ),
  main = indices[i],
  colpts = FinalDat$original.log[["Contam_col"]],
  boxwex = 0.5,
  ylab = ifelse(i == 1, "Richness (# of MOTUs)", ""),
  cex.axis = 0.9,
  ymin = 0,
  yscale = 6
)
}
dev.off()
