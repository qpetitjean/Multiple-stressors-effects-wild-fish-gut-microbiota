source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GutMicrobiome/BxpltFunc.R")

##############################################################################################################
#  Test homogeneity of variance among treatments group using betadisper                                      #
##############################################################################################################

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs"

# directory where the output files should be saved
DirSave <- file.path(savingDir, "Normalized_Data")

#######################
# Import the data     #
#######################

# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "Varia_contam_fullv5.csv"), dec = ".", sep = ";")

# import the distance matrix according to the beta diversity indices selected 
DistList <- readRDS(file.path(savingDir, "Normalized_Data/DistMatrices.RDS"))

# import the scaled coordinates of the samples in a 2d space (PCOA or NMDS), derived from the distance matrix
BetaDivCoords <- readRDS(file.path(savingDir, "Normalized_Data/BetaDivCoords.RDS"))

# import the results of the LMM approach on beta diversity indices
LMMRes <- readRDS(file.path(DirSave, "Signif_Effects", "LMMRes.rds"))

#########################################################################################
# compute betadisper and run permutation test over the normalized data                  #
#########################################################################################

# initialize an empty list to store the results of betadisper and permutation test
betadispRes <- list()

for(i in names(LMMRes)){
  for(k in names(LMMRes[[i]])){
    for(l in names(LMMRes[[i]][[k]])){
# merge the BetaDiv coords to the full dataset of the experiment
TempDist <- as.data.frame(BetaDivCoords[[l]][[k]][[i]][[ifelse(l == "PCOA", "vectors", "points")]])
TempDist[["Samples"]] <- rownames(TempDist)
TempDist <- merge(
  x = FullDat,
  y = TempDist,
  by.x = "Ind",
  by.y = "Samples",
  all = TRUE
)

# retrieve the significant effect from LMM results
toKeep <-
  unique(c(
    rownames(LMMRes[[i]][[k]][[l]][[ifelse(l == "PCOA", "Axis.1", "MDS1")]]$ModAnov)[-1],
    rownames(LMMRes[[i]][[k]][[l]][[ifelse(l == "PCOA", "Axis.2", "MDS2")]]$ModAnov)[-1]
  ))

# iterate trough significant effect highlighted in the LMM approach
for(j in toKeep) {
  if (is.character(TempDist[[j]]) | is.factor(TempDist[[j]])) {
    groups <-
      TempDist[match(labels(DistList[[i]][[k]]), TempDist[["Ind"]]), ]
    
    betadispRes[[i]][[k]][[j]][["betaDisper"]] <-
      vegan::betadisper(
        DistList[[i]][[k]],
        groups[[j]],
        type = "median",
        add = "lingoes",
        bias.adjust = TRUE
      )
    
    betadispRes[[i]][[k]][[j]][["PermTest"]] <-
      vegan::permutest(
        betadispRes[[i]][[k]][[j]][["betaDisper"]],
        pairwise = TRUE,
        permutations = 999,
        by = "margin"
      )
    
    betadispRes[[i]][[k]][[j]][["FDRadjustedP"]] <-
      p.adjust(betadispRes[[i]][[k]][[j]][["PermTest"]]$pairwise$observed, method = "fdr")
    
    TempDistOrd <- TempDist[order(TempDist[[j]]),]
    
    tempDf <-
      data.frame(
        gr = betadispRes[[i]][[k]][[j]][["betaDisper"]][["group"]],
        dist = betadispRes[[i]][[k]][[j]][["betaDisper"]][["distances"]],
        Ind = rownames(betadispRes[[i]][[k]][[j]][["betaDisper"]][["vectors"]])
      )
    tempDfPlot <- merge(
      x = tempDf,
      y = TempDistOrd,
      by.x = "Ind",
      by.y = "Ind",
      all = TRUE
    )
    
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots"))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots"))
    }
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots", i))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots", i))
    }
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots", i, k))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots", i, k))
    }
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots", i, k, l))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots", i, k, l))
    }
    ## save the figs representing distance to spatial median accroding to the treatments
    toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots", "BetadisperPlots", i, k, l, paste0(paste("Betadisper-Plots", i, k, l, j, sep = "_"), ".tif" ))
    if(file.exists(toSave)){
      unlink(toSave)
    } 
    if(j == "Contam"){
    newOrder <- c("NC", "C")
    } else if(j == "Inj") {
      newOrder <- c("PBS", "AMIX")
    } else {
      newOrder <- levels(factor(tempDfPlot[[j]]))
    } 
    tempDfPlot[["Ordered"]] <- factor(tempDfPlot[[j]], levels = newOrder)
    grDevices::tiff(filename = toSave,
                    width = 1200,
                    height = 1200,
                    res = 300,
                    compression = "lzw",
                    pointsize = 12)
    
    Bxplt(
      x = tempDfPlot[["Ordered"]],
      y =  tempDfPlot[["dist"]],
      fill = tempDfPlot[match(levels(tempDfPlot[["Ordered"]]), tempDfPlot[[j]]), paste(j, "col", sep = "_")],
      colpts = tempDfPlot[[paste(j, "col", sep = "_")]],
      las = ifelse(length(unique(tempDfPlot[[j]])) > 4, 2, 1),
      boxwex = 0.5,
      ylab = "Distance to spatial median",
      cex.axis = ifelse(length(unique(tempDfPlot[[j]])) > 4, 0.7, 1.1),
      main = paste("BetaDisper", paste(i, l, sep = " "), sep = " : "),
      yscale = 6)

    dev.off()
    
  } else if (is.numeric(TempDistOrd[[j]])) {
    next()
  }
}
    }
  }
}

saveRDS(betadispRes, file = file.path(DirSave, "Signif_Effects", "BetaDispRes.rds"), compress = TRUE)


################################################################
# summarize the results in html table for each normalization   #
################################################################

BestMods <- data.frame(
  mods = names(betadispRes),
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

# retrieve information from model summary and anova table
options(scipen = 999)
VRList <- list()
distances <- names(betadispRes[[1]])

for(Ord in c("PCOA", "NMDS")){
  Val <- setNames(lapply(distances, function(i) {
    ValuesRes <- data.frame()
    for (r in seq(nrow(BestMods))) {
      VarTest <- betadispRes[[BestMods[r, "mods"]]][[i]]
      effects <- list()
      if(is.null(names(VarTest))){
        next
      }
      for(v in names(VarTest)){
        effects[[v]] <- betadispRes[[BestMods[r, "mods"]]][[i]][[v]]$PermTest$tab
        rownames(effects[[v]])[1] <- v
      }
      Values <- do.call(rbind, effects)
      Values <- signif(Values, digits = 3)
      Values <- cbind(var = rownames(Values), Values)
      rownames(Values) <- NULL
      splitVar <- strsplit(Values$var, "\\.")
      Values["var"] <- sapply(splitVar, "[", 2)
      Values["mod"] <- sapply(splitVar, "[", 1)
      Values <- cbind(mods = rep(BestMods[r, "mods"], nrow(Values)), Values)
      ValuesRes <- rbind(ValuesRes, Values, stringsAsFactors = FALSE)
      rownames(ValuesRes) <- NULL
    }
    # leave the two first column without names
    names(ValuesRes)[c(1, 2)] <- c("", "")
    
    # specify the new predictor names (to display in the table)
    PredictorsNames <- data.frame(toreplace = unique(ValuesRes[[2]]), replacement = NA)
    PredictorsNames$replacement[PredictorsNames$toreplace == "Inj"] <- "Imm. Chall."
    PredictorsNames$replacement[PredictorsNames$toreplace == "Sexe"] <- "Sex"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Taille_mm_M"] <- "Size (cm)"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam"] <- "Contamination"
    PredictorsNames$replacement[PredictorsNames$toreplace == "ContPop"] <- "Origin"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:ContPop"] <- "Contam.:Origin"
    PredictorsNames$replacement[PredictorsNames$toreplace == "Contam:Inj"] <- "Contam.:Imm. Chall."
    PredictorsNames <- PredictorsNames[!is.na(PredictorsNames$replacement),]
    
    # rename the predictors
    for (pn in PredictorsNames[["toreplace"]]) {
      ValuesRes[[2]][which(ValuesRes[[2]] == pn)] <-
        PredictorsNames[["replacement"]][which(PredictorsNames[["toreplace"]] == pn)]
    }

    # replace p.value with large decimal by approximations
    ValuesRes[["Pr(>F)"]][which(ValuesRes[["Pr(>F)"]] < 0.0001)] <-
      "<0.0001"
    ValuesRes[["Pr(>F)"]][which(ValuesRes[["Pr(>F)"]] > 0.0001 &
                                   ValuesRes[["Pr(>F)"]] < 0.01)] <- "<0.01"
    ValuesRes[["Pr(>F)"]][which(ValuesRes[["Pr(>F)"]] > 0.01 &
                                   ValuesRes[["Pr(>F)"]] < 0.05)] <- "<0.05"

    return(ValuesRes)
  }), distances)
  VRList[[Ord]] <- Val
}
options(scipen = 0)

# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalisation methods)
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "BetaDispBestModTab"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "BetaDispBestModTab"))
} 


for(Ord in c("PCOA", "NMDS")){
  for(m in distances){
    Values <- VRList[[Ord]]
    ValIndex <- which(names(VRList[[Ord]]) == m)
    toTab <- Values[[ValIndex]]
    toTab[is.na(toTab)] <- ""
    l <- vector("list", nrow(unique(toTab[1])))
    l[[1]] <- kableExtra::kable_classic_2(
      kableExtra::kbl(toTab[-c(1, ncol(toTab))], 
                      align = "c", 
                      caption = paste0("<span style='font-size:12px; font-weight: bold; font-style: italic'>", 
                                       m, 
                                       " ",
                                       paste(Ord, sep = ": "), 
                                       "</span>")),
      bootstrap_options = c("striped", "hover", "condensed", "responsive"),
      full_width = F,
      html_font = "arial",
      font_size = 10
    )
    
    BestModsTemp <- BestMods[which(BestMods$mods %in% unique(toTab[[1]])),]
    for (i in seq(nrow(unique(toTab[1])))) {
      l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                          BestModsTemp[i, "modsId"],
                                          min(which(toTab[-ncol(toTab)][[1]] ==  BestModsTemp[i, "mods"])),
                                          max(which(toTab[-ncol(toTab)][[1]] ==  BestModsTemp[i, "mods"])),
                                          label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
    }
    l[[length(l)]]
    kableExtra::save_kable(l[[length(l)]], 
                           file.path(DirSave, "Signif_Effects", "BetaDispBestModTab", 
                                     paste0(paste("BetaDispBestModTab", m, Ord, sep = "_"), 
                                            ".HTML")))
  }
}

# generate the table including all metric for the original data log transformed (Table ...)
## retrieve the summary from the distances of the original Log transformed data in PCOA ordination (for both axes)
dataNorm <- "original.log"
ordination <- "PCOA"

ValuesPart <- VRList[[ordination]]
OrLog <- lapply(ValuesPart, function(x)
  x[which(x[[1]] == dataNorm), ])
OrLogDf <- do.call(rbind,OrLog)
OrLogDf$index <- gsub("\\.\\d+", "", rownames(OrLogDf))
rownames(OrLogDf) <- NULL
OrLogDf[is.na(OrLogDf)] <- ""

l <- vector("list", length(unique(OrLogDf[["index"]])))
l[[1]] <- kableExtra::kable_classic_2(
  kableExtra::kbl(
    OrLogDf[-c(1, length(OrLogDf),length(OrLogDf)-1 )],
    align = "c",
    caption = paste0(
      "<span style='font-size:12px; font-weight: bold; font-style: italic'>",
      paste(dataNorm,
            ordination, sep = " "),
      "</span>"
    )
  ),
  bootstrap_options = c("striped", "hover", "condensed", "responsive"),
  full_width = F,
  html_font = "arial",
  font_size = 10
)

for (i in seq_along(distances)) {
  l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                      distances[i], 
                                      min(which(OrLogDf[["index"]] == distances[i])),
                                      max(which(OrLogDf[["index"]] == distances[i])),
                                      label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
}
kableExtra::save_kable(l[[length(l)]], 
                       file.path(DirSave, "Signif_Effects", "BetaDispBestModTab", 
                                 paste0(paste("BetaDispBestModTab", dataNorm, ordination, sep = "_"), 
                                        ".HTML")))

# compute mean distance to median and se for some indices that give significant differences 
## For "Robust.Aitchison" in contamination treatment
df <- data.frame(dist = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Contam$betaDisper$distances,
                 group = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Contam$betaDisper$group)

meanRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) mean(x, na.rm = T))
seRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) sd(x, na.rm = T))$x /
  stats::aggregate(df$dist, list(df$group), function(x)
    sqrt(sum(!is.na(x))))$x

## For "Robust.Aitchison" in contamination treatment
df <- data.frame(dist = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Sexe$betaDisper$distances,
                 group = betadispRes[[dataNorm]][["Robust.Aitchison"]]$Sexe$betaDisper$group)

meanRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) mean(x, na.rm = T))
seRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) sd(x, na.rm = T))$x /
  stats::aggregate(df$dist, list(df$group), function(x)
    sqrt(sum(!is.na(x))))$x

p <- betadispRes[[dataNorm]][["Robust.Aitchison"]]$Sexe$FDRadjustedP

## For "Hellinger" in contamination treatment
df <- data.frame(dist = betadispRes[[dataNorm]][["Hellinger"]]$Sexe$betaDisper$distances,
                 group = betadispRes[[dataNorm]][["Hellinger"]]$Sexe$betaDisper$group)

meanRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) mean(x, na.rm = T))
seRA <- stats::aggregate(df$dist, list(df$group), FUN = function(x) sd(x, na.rm = T))$x /
  stats::aggregate(df$dist, list(df$group), function(x)
    sqrt(sum(!is.na(x))))$x

p <- betadispRes[[dataNorm]][["Hellinger"]]$Sexe$FDRadjustedP
