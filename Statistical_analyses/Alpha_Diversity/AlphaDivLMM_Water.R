##############################################################################################################
#  Compute Alpha diversity indices and test treatments effects using LMM and dredge method                   #
#############################################################################################################
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/multiNorm.R") # import a function used to compute multiple normalization methods
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/DredgedLMM.R")  # import a function used to retrieve the summary, anova table and Rsquared of the most straightforward mixed effect model using rank selection (e.g., AICc)
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/alphaDiv.R") # import a function compute alpha diversity indices from a matrix containing read counts per sample
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GutMicrobiome/BxpltFunc.R") # import a function compute draw customized boxplots

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported and saved
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs"

#######################
# Import the data     #
#######################

# import the phylogenetic tree (.nwk)
PhyloTree <-
  ape::read.tree(file.path(savingDir, "PhyloTree.nwk"))


# import the cleaned dataset (metabaR object) and select only gut samples
labWater <- readRDS(file.path(savingDir, "Preprocessing-Metabar/fguts_Bact_agg_MergedRep.RDS"))
labWater <- metabaR::subset_metabarlist(labWater,
                                       table = "pcrs",
                                       indices = labWater$pcrs$matrix == "water")

names(labWater$samples)[which(names(labWater$samples) == "treatment_contam")] <- "Contam"
names(labWater$samples)[which(names(labWater$samples) == "pop")] <- "Pop"
names(labWater$samples)[which(names(labWater$samples) == "treatment_inj")] <- "Inj"

# add color for treatment groups 
## manually create interaction treatment and specify colors to ease graphical representation
labWater$samples[["Contam:Pop"]] <-
  paste(labWater$samples[["Contam"]], labWater$samples[["Pop"]], sep = "_")
labWater$samples[["Contam:Inj"]] <-
  paste(labWater$samples[["Contam"]], labWater$samples[["Inj"]], sep = "_")
labWater$samples[["Inj:Pop"]] <-
  paste(labWater$samples[["Inj"]], labWater$samples[["Pop"]], sep = "_")

labWater$samples[["Contam_col"]] <- NA
labWater$samples[["Contam_col"]][grep("NC", labWater$samples[["Contam"]])] <- adjustcolor("#006400", alpha.f = 0.3)
labWater$samples[["Contam_col"]][grep("^((?!NC).)*$", labWater$samples[["Contam"]], perl = T)] <- adjustcolor("#B22222", alpha.f = 0.3)

labWater$samples[["Inj_col"]] <- NA
labWater$samples[["Inj_col"]][grep("PBS", labWater$samples[["Inj"]])] <- adjustcolor("#808080", alpha.f = 0.3)
labWater$samples[["Inj_col"]][grep("LPS", labWater$samples[["Inj"]], perl = T)] <- adjustcolor("#470C8E", alpha.f = 0.3)

labWater$samples[["Pop_col"]] <- NA
labWater$samples[["Pop_col"]][grep("ELEV", labWater$samples[["Pop"]])] <- adjustcolor("#09475e", alpha.f = 0.3)
labWater$samples[["Pop_col"]][grep("ARIMAS", labWater$samples[["Pop"]], perl = T)] <- adjustcolor("#40457f", alpha.f = 0.3)
labWater$samples[["Pop_col"]][grep("CELCAB", labWater$samples[["Pop"]], perl = T)] <- adjustcolor("#B0C4DE", alpha.f = 0.3)
labWater$samples[["Pop_col"]][grep("AUSCOR", labWater$samples[["Pop"]], perl = T)] <- adjustcolor("#91296b", alpha.f = 0.3)
labWater$samples[["Pop_col"]][grep("RIOU", labWater$samples[["Pop"]], perl = T)] <- adjustcolor("#b22222", alpha.f = 0.3)

labWater$samples[["Contam:Pop_col"]] <- NA
labWater$samples[["Contam:Pop_col"]][grep("NC", labWater$samples[["Contam"]])] <- adjustcolor("#006400", alpha.f = 0.3)
labWater$samples[["Contam:Pop_col"]][grep("^((?!NC).)*$", labWater$samples[["Contam"]], perl = T)] <- adjustcolor("#B22222", alpha.f = 0.3)

labWater$samples[["Contam:Inj_col"]] <- NA
labWater$samples[["Contam:Inj_col"]][grep("NC", labWater$samples[["Contam"]])] <- adjustcolor("#006400", alpha.f = 0.3)
labWater$samples[["Contam:Inj_col"]][grep("^((?!NC).)*$", labWater$samples[["Contam"]], perl = T)] <- adjustcolor("#B22222", alpha.f = 0.3)


#######################################
# compute multiple normalization      #
#######################################
normalizedDat <- multiNorm(
  dataList = labWater,
  norm = c(
    "original",
    "original.log",
    "read.prop",
    "read.prop.log",
    "read.CSS",
    "read.CSS.log",
    "read.rare",
    "read.rare.log",
   #"read.deseq",
   #"read.deseq.log",
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
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "LMM"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water","LMM"))
} 
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots","Water", "LMM", "AlphaDivPlots"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots","Water", "LMM", "AlphaDivPlots"))
} 


# create several versions of the full dataset by merging the alphadiv indices to the full dataset of the experiment
FinalDat <- lapply(alphaDivindices, function(x)
  merge(
    x = labWater$samples,
    y = x,
    by.x = "sample_id",
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
   
    ModL <- lm(Temp[[indices[[j]]]] ~ Contam + Pop + Inj, data = Temp)
    ModSum <- summary(ModL)
    ModAnov <- car::Anova(ModL, type = "3")
    res[[names(FinalDat)[i]]][[indices[j]]] <- list(ModL = ModL, ModSum = ModSum, ModAnov = ModAnov)
    # save performance plot (LMM diagnostic)
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "LMM", "AlphaDivPlots", names(FinalDat)[i]))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "LMM", "AlphaDivPlots", names(FinalDat)[i]))
    } 
    toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots", "Water", "LMM", "AlphaDivPlots", names(FinalDat)[i], paste0(paste("Mod_Perf_AlphaDivLMM", names(FinalDat)[i], indices[j], sep = "_"), ".tif" ))
    if(file.exists(toSave)){
      unlink(toSave)
    }
    grDevices::tiff(filename = toSave,
                    width = 1200,
                    height = 1200,
                    res = 100,
                    compression = "lzw",
                    pointsize = 8)
    #print(performance::check_model(res[[names(FinalDat)[i]]][[indices[j]]]$ModL))
    hist(resid(res[[names(FinalDat)[i]]][[indices[j]]]$ModL))
    dev.off()
  }
}

# save statistical test outputs

if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Water"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Water"))
} 
saveRDS(res, file = file.path(DirSave, "Signif_Effects", "Water", "LMAlphaDivResWater.rds"))

## plot the results
for (i in seq_along(FinalDat)){
  Temp <- FinalDat[[i]]
  for(j in seq_along(indices)){
    if(length(which(is.na(Temp[indices[[j]]]))) > 0 && 
       nrow(Temp[-which(is.na(Temp[indices[[j]]])), ]) == 0){
      next
    }
    toSave <- file.path(DirSave, "Signif_Effects/Signif_Plots/Water/LMM/AlphaDivPlots", names(FinalDat)[i], gsub(":", "X", paste0(paste("AlphaDiv", names(FinalDat)[i], gsub("\\.", "_", indices[[j]]), sep = "_"), ".tif")))
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

##### plot differences among Populations
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

res <- readRDS(file.path(DirSave, "Signif_Effects", "Water", "LMAlphaDivResWater.rds"))
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
    #"DESEQ",
    #"DESEQ (log-transformed)",
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
        paste0("n = ", nrow(model.frame(res[[BestMods[r, "mods"]]][[i]]$ModL)))
      # Rsquared
      BestMods[r, paste("R2m", i, sep = "_")] <-
        paste0("Mult. R2 = ", signif(res[[BestMods[r, "mods"]]][[i]]$ModSum$r.squared, digits = 3))
      BestMods[r, paste("R2c", i, sep = "_")] <-
        paste0("Adj. R2 = ", signif(res[[BestMods[r, "mods"]]][[i]]$ModSum$adj.r.squared, digits = 3))
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
  SSq <- res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Sum Sq"]]
  F_value = res[[BestMods[r, "mods"]]][[i]]$ModAnov[["F value"]]
  p.value = res[[BestMods[r, "mods"]]][[i]]$ModAnov[["Pr(>F)"]]
  if(nrow(coef) != length(df)){
    match <- unlist(find_matches(rownames(coef)))
    # deal with interactions
    isCont <- grep("Cont", match)
    if(length(isCont)>0 & length(grep(":", isCont)) == 0 ){
      match <- match[-isCont]
    }
    for(m in seq(length(match)-1)){
    char <- substr(match[m], 1, 3)
    index <- grep(match[m], rownames(coef))
    df <- append(df, df[index], after=index)
    SSq <- append(SSq, SSq[index], after=index)
    F_value <- append(F_value, F_value[index], after=index)
    p.value <- append(p.value, p.value[index], after=index)
    }
  }
  coef <- coef[,-which(colnames(coef) == "Pr(>|t|)")]
  coef <- rbind(coef, Residuals = NA)
  Values <- cbind(coef, df, SSq, F_value, p.value)
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

Val
# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalization methods)
## create a directory to store the table containing LMM results for alpha diversity indices according to selected normalization
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "LMMAlphaDivWaterBestModTab"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "LMMAlphaDivWaterBestModTab"))
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
kableExtra::save_kable(l[[length(l)]], file.path(DirSave, "Signif_Effects", "LMMAlphaDivWaterBestModTab", paste0(paste("LMMAlphaDivWaterBestMod", m, sep = "_"), ".HTML")))
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
kableExtra::save_kable(l[[length(l)]], file.path(DirSave, "Signif_Effects", "LMMAlphaDivWaterBestModTab", paste0(
  paste("LMMAlphaDivWaterBestMod", dataNorm, sep = "_"), ".HTML"
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
toSave <- file.path(DirSave, "Signif_Effects/Signif_Plots/Water/LMM/AlphaDivPlots/BP_RICH_LogORDat.tif")
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
  yscale = 8
)
dev.off()


# display the boxplot of differences among Contamination treatment for the original log transformed data for all indices
# Reorder the factor levels
newOrder <- c("NC", "C")
FinalDat$original.log[["ContamOrd"]] <- factor(FinalDat$original.log[["Contam"]], levels = newOrder)
toSave <- file.path(DirSave, "Signif_Effects/Signif_Plots/Water/LMM/AlphaDivPlots/BP_AllIndices_LogORDat.tif")
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
