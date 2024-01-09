##############################################################################################################
#  Compute Beta diversity indices and test treatments effects using LMM and dredge method                   #
#############################################################################################################
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/multiNorm.R") # import a function used to compute multiple normalization methods
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/DredgedLMM.R")  # import a function used to retrieve the summary, anova table and Rsquared of the most straightforward mixed effect model using rank selection (e.g., AICc)
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/alphaDiv.R") # import a function compute alpha diversity indices from a matrix containing read counts per sample
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/multiDist.R") # import a function used to compute multiple distances
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GutMicrobiome/BxpltFunc.R") # import a function compute draw customized boxplots

#######################
# specify some path   #
#######################
# the path to the directory from where file should be both imported 
savingDir <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs"

# directory where the output files should be saved
DirSave <- "D:/POSTDOC_INP_GOLFECH_2023/Outputs/Normalized_Data"

#######################
# Import the data     #
#######################
# import the full dataset from the experiment (treatment & behavioral & physiological measures)
FullDat <-
  read.csv2(file.path(savingDir, "Varia_contam_fullv5.csv"), dec = ".", sep = ";")

# import the phylogenetic tree (.nwk)
PhyloTree <-
  ape::read.tree(file.path(savingDir, "PhyloTree.nwk"))

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

###############################################################
# import Picrust2 results
###############################################################
## import EC data including description column (tibble)
ECPath <-  as.data.frame(readr::read_tsv(file.path(savingDir, "Picrust2Output/EC_pathways_out/pred_metagenome_unstrat_descrip.tsv.gz")))

# manual change of the PWY-5182 and ARGORNPROST-PWY pathways which are not found in the Metacyc db 
# but found trough manual search using description given by picrust2 
ECPath[grep("PWY-5182", ECPath[,1]), 1] <- "TOLUENE-DEG-3-OH-PWY"
ECPath[grep("ARGORNPROST-PWY", ECPath[,1]), 1] <- "ARGDEGRAD-PWY"

ECMat <- ECPath[,-c(1,2)]
rownames(ECMat) <- ECPath[,1]
colnames(ECMat) <- labFish[["samples"]][["Ind"]][match(colnames(ECMat), rownames(labFish[["samples"]]))]

ECDat <- list(reads = ECMat, samples = labFish[["samples"]])

#######################################
# compute multiple normalization      #
#######################################
normalizedDat <- multiNorm(
  dataList = ECDat,
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
  ),
  TaxRow = T
)

#############################################################
# compute Beta div. distances over the normalized data      #
#############################################################

# the distances to compute 
distances <- c("Bray","Jaccard","Robust.Aitchison","Hellinger")

## Visualize Bray Curtis distances (PCOA)
# initialize some list to retrieve results and plots
DistList <- list()
BetaDivCoords <- list()
PlotListPop <- list()
PlotListContam <- list()
PlotListInj <- list()
LMMRes <-list()
EnvFitRes <- list()

# initialize progress bar
total = length(normalizedDat)
pb <-
  progress::progress_bar$new(format = "Processing [:bar] :current/:total (:percent)", total = total)
pb$tick(0)
for (i in names(normalizedDat)) {
  temp <- normalizedDat[[i]]
  
  # compute the specified distances 
  DistList[[i]] <- multiDist(as.matrix(temp), dist = distances, PhyloTree)
  
  # iterate trough the list of distance matrix to retrieve scaled coordinates of the samples in a 2d space (PCOA or NMDS)
  for (k in distances) {
    ## compute PCOA coords
    if (class(DistList[[i]][[k]]) == "dist") {
      BetaDivCoords[["PCOA"]][[k]][[i]] <-
        ape::pcoa(DistList[[i]][[k]], correction = "lingoes")
      ## compute NMDS coords
      BetaDivCoords[["NMDS"]][[k]][[i]] <-
        vegan::monoMDS(DistList[[i]][[k]], k = 2, model = "global")
    } else if (class(DistList[[i]][[k]]) == "poidist") {
      BetaDivCoords[["PCOA"]][[k]][[i]] <-
        ape::pcoa(DistList[[i]][[k]]$dd, correction = "lingoes")
      ## compute NMDS coords
      BetaDivCoords[["NMDS"]][[k]][[i]] <-
        vegan::monoMDS(DistList[[i]][[k]]$dd, k = 2, model = "global")
    } else{
      stop("class",
           class(DistList[[i]][[k]]),
           "not supported by ape::pcoa and vegan::monoMDS")
    }
    
    # retrieve the coordinates in PCOA and NMDS space and test the correlation between individual's scores and treatments and covariates. 
    # Also, test the significance of the treatments (AICc model selection) on the principal axis: axis 1 and 2
    for (l in c("PCOA", "NMDS")) {
      ### merge the BetaDiv coords to the full dataset of the experiment
      TempDist <- as.data.frame(BetaDivCoords[[l]][[k]][[i]][[ifelse(l == "PCOA", "vectors", "points")]])
      TempDist[["Samples"]] <- rownames(TempDist)
      TempDist <- merge(
        x = FullDat,
        y = TempDist,
        by.x = "Ind",
        by.y = "Samples",
        all = TRUE
      )
      
      ## run envfit analysis (check whether PCOA or NMDS scores are correlated with treatments and covariates)
      TempDistnoNA <- TempDist[!is.na(ifelse(rep(l == "PCOA", nrow(TempDist)), TempDist[["Axis.1"]], TempDist[["MDS1"]])),]
      axisLoc <- grep(ifelse(l == "PCOA", "Axis.", "MDS"), names(TempDistnoNA),  fixed = T)
      ord <- TempDistnoNA[,axisLoc]
      EnvFitRes[[i]][[k]][[l]] <- vegan::envfit(ord ~ Contam + Inj + ContPop + Sexe + Taille_mm_M, TempDistnoNA, permutations = 999)
      
      for(axe in ifelse(rep(l == "PCOA", 2), c("Axis.1", "Axis.2"), c("MDS1", "MDS2"))){        
        TempDist[["Bac"]] <- as.character(TempDist[["Bac"]])
        ## dredge the full model and retrieve the summary, anova table and R2 of the most straightforward one 
        BestMod <- DredgedLMM(TempDist, 
                              expVar = axe,
                              respExpr = "(Contam + ContPop + Inj)^3 + Taille_mm_M + Sexe",
                              random = "(1 |Session / Bac) + (1 | Pop)",
                              method = "weight",
                              rank = "AICc")
        ## store the result in a list
        LMMRes[[i]][[k]][[l]][[axe]] <- BestMod
        ##### create output folder to save figs
        if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions"))) == 0) {
          dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions"))
        }
        if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM"))) == 0) {
          dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM"))
        }
        if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots"))) == 0) {
          dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots"))
        }
        if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots", i))) == 0) {
          dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots", i))
        }
        if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots", i, k))) == 0) {
          dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots", i, k))
        }
        if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots", i, k, l))) == 0) {
          dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots", i, k, l))
        }
        toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots", "Functions", "LMM", "BetaDivPlots", i, k, l, paste0("Mod_Perf_", axe, ".tif" ))
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
        
        toKeep <- rownames(BestMod$ModAnov)[-1]
        
        if(length(toKeep) == 0){
          next
        }
        ##### save the figs representing difference of PCOA Scores (composition among significant treatments)
        toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots","Functions", "LMM", "BetaDivPlots", i, k, l, paste0("Signif-Effects-Plots_", axe, ".tif" ))
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
          TempDistOrd <- TempDist[order(TempDist[[keepVar]]), ]
          if(is.character(TempDistOrd[[keepVar]]) | is.factor(TempDistOrd[[keepVar]])){
            Bxplt(
              x = TempDistOrd[[keepVar]],
              y =  TempDistOrd[[axe]],
              fill = TempDistOrd[match(unique(TempDistOrd[[keepVar]]), TempDistOrd[[keepVar]]), paste(keepVar, "col", sep = "_")],
              colpts = TempDistOrd[[paste(keepVar, "col", sep = "_")]],
              las = ifelse(length(unique(TempDistOrd[[keepVar]])) > 4, 2, 1),
              boxwex = 0.5,
              ylab = paste(l, axe, sep = " - "),
              cex.axis = ifelse(length(unique(TempDistOrd[[keepVar]])) > 4, 0.7, 1.1),
              main = paste(i, keepVar, sep = " : "),
              yscale = 6)
          }else if(is.numeric(TempDistOrd[[j]])) {
            plot(x = TempDistOrd[[keepVar]],
                 y =  TempDistOrd[[axe]],
                 xlab = keepVar,
                 ylab = paste(l, axe, sep = " - "),
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
      
      ### store the plot (PCOA or NMDS)
      if(l == "PCOA") {
        VarAxis1 <- (BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]][[1]]/
                       sum(BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]])) * 100
        VarAxis2 <- (BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]][[2]]/
                       sum(BetaDivCoords[[l]][[k]][[i]][["values"]][["Eigenvalues"]])) * 100
      }
      
      ordLab <- c("ELEV", "ARIMAS", "CELCAB", "AUSCOR", "RIOU")
      ordCol <- as.character(TempDist$Pop_col[match(ordLab, TempDist$Pop)])
      PlotListPop[[l]][[k]][[i]] <-
        ggplot2::ggplot(TempDist,
                        ggplot2::aes(x = .data[[ifelse(l == "PCOA", "Axis.1", "MDS1")]], 
                                     y = .data[[ifelse(l == "PCOA", "Axis.2", "MDS2")]], 
                                     colour = .data[["Pop_col"]], 
                                     fill = .data[["Pop_col"]])
        ) +
        ggplot2::geom_point(size = 4) +
        ggplot2::ggtitle(ifelse(l == "PCOA", i, paste0(
          i,
          " (stress = ",
          signif(BetaDivCoords[[l]][[k]][[i]][["stress"]], digits = 2),
          ")"
        ))) +
        ggplot2::xlab(ifelse(l == "PCOA", paste0("Axis.1 (", round(VarAxis1, digits = 1), "%", ")"), "MDS1")) +
        ggplot2::ylab(ifelse(l == "PCOA", paste0("Axis.2 (", round(VarAxis2, digits = 1), "%", ")"), "MDS2")) +
        ggplot2::stat_ellipse(geom = "polygon",
                              ggplot2::aes(fill = .data[["Pop_col"]], colour = .data[["Pop_col"]]), 
                              alpha = 0.25,
                              type = "t") +
        ggplot2::scale_colour_identity(guide = "legend", 
                                       breaks = ordCol, 
                                       labels = ordLab, 
                                       name = "Population") + 
        ggplot2::scale_fill_identity(guide = "legend", 
                                     breaks = ordCol, 
                                     labels = ordLab, 
                                     name = "Population")+
        ggplot2::theme_classic()
      
      ordLab <- c("NC", "C")
      ordCol <- as.character(TempDist$Contam_col[match(ordLab, TempDist$Contam)])
      PlotListContam[[l]][[k]][[i]] <-
        ggplot2::ggplot(TempDist,
                        ggplot2::aes(x = .data[[ifelse(l == "PCOA", "Axis.1", "MDS1")]], 
                                     y = .data[[ifelse(l == "PCOA", "Axis.2", "MDS2")]], 
                                     colour = .data[["Contam_col"]], 
                                     fill = .data[["Contam_col"]])
        ) +
        ggplot2::geom_point(size = 4) +
        ggplot2::ggtitle(ifelse(l == "PCOA", i, paste0(
          i,
          " (stress = ",
          signif(BetaDivCoords[[l]][[k]][[i]][["stress"]], digits = 2),
          ")"
        ))) +
        ggplot2::xlab(ifelse(l == "PCOA", paste0("Axis.1 (", round(VarAxis1, digits = 1), "%", ")"), "MDS1")) +
        ggplot2::ylab(ifelse(l == "PCOA", paste0("Axis.2 (", round(VarAxis2, digits = 1), "%", ")"), "MDS2")) +
        ggplot2::stat_ellipse(geom = "polygon",
                              ggplot2::aes(fill = .data[["Contam_col"]], colour = .data[["Contam_col"]]), 
                              alpha = 0.25,
                              type = "t") + 
        ggplot2::scale_colour_identity(guide = "legend", 
                                       breaks = ordCol, 
                                       labels = ordLab, 
                                       name = "Contamination") + 
        ggplot2::scale_fill_identity(guide = "legend", 
                                     breaks = ordCol, 
                                     labels = ordLab, 
                                     name = "Contamination")+
        ggplot2::theme_classic()
      
      ordLab <- c("PBS", "AMIX")
      ordCol <- as.character(TempDist$Inj_col[match(ordLab, TempDist$Inj)])
      PlotListInj[[l]][[k]][[i]] <-
        ggplot2::ggplot(TempDist,
                        ggplot2::aes(x = .data[[ifelse(l == "PCOA", "Axis.1", "MDS1")]], 
                                     y = .data[[ifelse(l == "PCOA", "Axis.2", "MDS2")]], 
                                     colour = .data[["Inj_col"]], 
                                     fill = .data[["Inj_col"]])
        ) +
        ggplot2::geom_point(size = 4) +
        ggplot2::ggtitle(ifelse(l == "PCOA", i, paste0(
          i,
          " (stress = ",
          signif(BetaDivCoords[[l]][[k]][[i]][["stress"]], digits = 2),
          ")"
        ))) +
        ggplot2::xlab(ifelse(l == "PCOA", paste0("Axis.1 (", round(VarAxis1, digits = 1), "%", ")"), "MDS1")) +
        ggplot2::ylab(ifelse(l == "PCOA", paste0("Axis.2 (", round(VarAxis2, digits = 1), "%", ")"), "MDS2")) +
        ggplot2::stat_ellipse(geom = "polygon",
                              ggplot2::aes(fill = .data[["Inj_col"]], colour = .data[["Inj_col"]]), 
                              alpha = 0.25,
                              type = "t") +
        ggplot2::scale_colour_identity(guide = "legend", 
                                       breaks = ordCol, 
                                       labels = ordLab, 
                                       name = "Imm. Chall.") + 
        ggplot2::scale_fill_identity(guide = "legend", 
                                     breaks = ordCol, 
                                     labels = ordLab, 
                                     name = "Imm. Chall.")+
        ggplot2::theme_classic()
    }
  }
  pb$tick(1)
}

for(h in c("Pop", "Contam", "Inj")) {
  PlotList <- get(paste0("PlotList", h), envir = .GlobalEnv)
  for (i in names(PlotList)) {
    for (j in names(PlotList[[i]])) {
      PlotsToSave <-
        ggpubr::ggarrange(plotlist = PlotList[[i]][[j]],
                          ncol = 4,
                          nrow = 3)
      if (length(list.dirs(file.path(DirSave, "Functions", i))) == 0) {
        dir.create(file.path(DirSave, "Functions", i))
      }
      if (length(list.dirs(file.path(DirSave, "Functions", i, j))) == 0) {
        dir.create(file.path(DirSave, "Functions", i, j))
      }
      toSave <- paste0(file.path(
        DirSave, "Functions", i, j, paste(h, i, j, "Functions", "Norm", sep = "_")
      ), ".svg")
      if(file.exists(toSave)){
        unlink(toSave)
      } 
      ggplot2::ggsave(
        file = toSave,
        plot = PlotsToSave,
        width = 16,
        height = 12
      )
    }
  }
}

# save the distances matrices as .RDS
saveRDS(DistList, paste0(file.path(DirSave, "DistMatrices_Functions"), ".RDS"), compress = TRUE)

# save the BetaDivCoords (PCOA and NMDS) as .RDS
saveRDS(BetaDivCoords, paste0(file.path(DirSave, "BetaDivCoords_Functions"), ".RDS"), compress = TRUE)

# save statistical test outputs
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Functions"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "Functions"))
}
saveRDS(LMMRes, file = file.path(DirSave, "Signif_Effects", "Functions", "LMMRes_Functions.rds"), compress = TRUE)

# save Envfit test outputs
saveRDS(EnvFitRes, file = file.path(DirSave, "Signif_Effects", "Functions", "EnvfitRes_Functions.rds"), compress = TRUE)

################################################################
# summarize the results in html table for each normalization   #
################################################################

##################
# For LMM        #
#################

# make the statistical summary table (PCOA axis 1 and 2 LMM + EnvFit Res for all distances) with original log transformed data 
## retrieve the sample size and the Rsquared for each ordination and axis
BMList <- list()
for(Ord in c("PCOA", "NMDS")){
  for(axe in ifelse(rep(Ord == "PCOA", 2), c("Axis.1", "Axis.2"), c("MDS1", "MDS2"))){
    # create a temporary table grouping the names of the best models and the name of the corresponding response variable
    BestMods <- data.frame(
      mods = names(LMMRes),
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
    for(i in distances){
      BestMods[[paste("n", i, sep = "_")]] <- rep(NA, nrow(BestMods))
      BestMods[[paste("R2m", i, sep = "_")]] <- rep(NA, nrow(BestMods))
      BestMods[[paste("R2c", i, sep = "_")]] <- rep(NA, nrow(BestMods))
      for (r in seq(nrow(BestMods))) {
        # n
        BestMods[r, paste("n", i, sep = "_")] <-
          paste0("n = ", LMMRes[[BestMods[r, "mods"]]][[i]][[Ord]][[axe]]$ModSum$devcomp$dims[["n"]])
        # Rsquared
        BestMods[r, paste("R2m", i, sep = "_")] <-
          paste0("R2m = ", signif(LMMRes[[BestMods[r, "mods"]]][[i]][[Ord]][[axe]]$Rsquared[[1]], digits = 3))
        BestMods[r, paste("R2c", i, sep = "_")] <-
          paste0("R2c = ", signif(LMMRes[[BestMods[r, "mods"]]][[i]][[Ord]][[axe]]$Rsquared[[2]], digits = 3))
      }
    }
    BMList[[Ord]][[axe]] <- BestMods
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
VRList <- list()

for(Ord in c("PCOA", "NMDS")){
  for(axe in ifelse(rep(Ord == "PCOA", 2), c("Axis.1", "Axis.2"), c("MDS1", "MDS2"))){
    Val <- setNames(lapply(distances, function(i) {
      ValuesRes <- data.frame()
      for (r in seq(nrow(BestMods))) {
        coef <- LMMRes[[BestMods[r, "mods"]]][[i]][[Ord]][[axe]]$ModSum[["coefficients"]]
        df <- LMMRes[[BestMods[r, "mods"]]][[i]][[Ord]][[axe]]$ModAnov[["Df"]]
        Chisq = LMMRes[[BestMods[r, "mods"]]][[i]][[Ord]][[axe]]$ModAnov[["Chisq"]]
        p.value = LMMRes[[BestMods[r, "mods"]]][[i]][[Ord]][[axe]]$ModAnov[["Pr(>Chisq)"]]
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
    }), distances)
    VRList[[Ord]][[axe]] <- Val
  }
}
options(scipen = 0)

# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalisation methods)
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "LMMBetaDivFunctionsBestModTab"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "LMMBetaDivFunctionsBestModTab"))
} 

for(Ord in c("PCOA", "NMDS")){
  for(axe in ifelse(rep(Ord == "PCOA", 2), c("Axis.1", "Axis.2"), c("MDS1", "MDS2"))){
    for(m in distances){
      BestMods <- BMList[[Ord]][[axe]]
      Values <- VRList[[Ord]][[axe]]
      ValIndex <- which(names(VRList[[Ord]][[axe]]) == m)
      l <- vector("list", nrow(BestMods))
      l[[1]] <- kableExtra::kable_classic_2(
        kableExtra::kbl(Values[[ValIndex]][-1], 
                        align = "c", 
                        caption = paste0("<span style='font-size:12px; font-weight: bold; font-style: italic'>", 
                                         m, 
                                         " ",
                                         paste(Ord, axe, sep = ": "), 
                                         "</span>")),
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
                                            min(which(Values[[ValIndex]][[1]] ==  BestMods[i, "mods"])),
                                            max(which(Values[[ValIndex]][[1]] ==  BestMods[i, "mods"])),
                                            label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
      }
      l[[length(l)]]
      kableExtra::save_kable(l[[length(l)]], 
                             file.path(DirSave, "Signif_Effects", "LMMBetaDivFunctionsBestModTab", 
                                       paste0(paste("LMMBetaDivFunctionsBestModTab", m, Ord, axe, sep = "_"), 
                                              ".HTML")))
    }
  }
}
# generate the table including all metric for the original data log transformed (Table 3 & 4)
## retrieve the summary from the distances of the original Log transformed data in PCOA ordination (for both axes)
dataNorm <- "original.log"
ordination <- "PCOA"

TabList <- list()

for(axis in c("Axis.1", "Axis.2")){
  BestModPart <- BMList[[ordination]][[axis]]
  BestModPart <- BestModPart[which(BestModPart$mods == dataNorm),]
  BestModPart <- BestModPart[,-c(1,2)]
  BestModPart <- setNames(split.default(BestModPart, gl(ncol(BestModPart) / 3, 3)), distances)
  ValuesPart <- VRList[[ordination]][[axis]]
  
  OrLog <- lapply(ValuesPart, function(x)
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
        paste(dataNorm,
              ordination,
              axis, sep = " "),
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
  TabList[[axis]] <- l[[length(l)]]
  kableExtra::save_kable(TabList[[axis]], 
                         file.path(DirSave, "Signif_Effects", "LMMBetaDivFunctionsBestModTab", 
                                   paste0(paste("LMMBetaDivFunctionsBestModTab", dataNorm, ordination, axis, sep = "_"), 
                                          ".HTML")))
}


# save the corresponding PCOA for contamination treatment (the only one significant treatment) (Figure 1B)
PlotListContam[[ordination]][["Robust.Aitchison"]][[dataNorm]]
if (length(list.dirs(file.path(DirSave, "Functions"))) == 0) {
  dir.create(file.path(DirSave,"Functions"))
} 
toSave <- paste0(file.path(
  DirSave, "Functions", ordination, "Robust.Aitchison", paste("Contam", ordination, dataNorm, "Robust.Aitchison", sep = "_")
), ".svg")
if(file.exists(toSave)){
  unlink(toSave)
} 

# change the size of the dots and remove the title
PlotListContam[[ordination]][["Robust.Aitchison"]][[dataNorm]]$layers[[1]]$aes_params$size <- 2
PlotListContam[[ordination]][["Robust.Aitchison"]][[dataNorm]]$labels$title <- ""

ggplot2::ggsave(
  file = toSave,
  plot = PlotListContam[[ordination]][["Robust.Aitchison"]][[dataNorm]],
  dpi=200,
  width = 5,
  height = 3
)


##################
# For EnvFit     #
#################

## make and save the table for the Envfit Results (one table per distance and per ordination e.g. Bray PCOA)
# create a list of normalisation and rename them
BestMods <- data.frame(
  mods = names(EnvFitRes),
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

for(Ord in c("PCOA", "NMDS")){
  Val <- setNames(lapply(distances, function(i) {
    ValuesRes <- data.frame()
    for (r in seq(nrow(BestMods))) {
      RVec <- EnvFitRes[[BestMods[r, "mods"]]][[i]][[Ord]]$vectors[["r"]]
      pVec <- EnvFitRes[[BestMods[r, "mods"]]][[i]][[Ord]]$vectors[["pvals"]]
      RFac <- EnvFitRes[[BestMods[r, "mods"]]][[i]][[Ord]]$factors[["r"]]
      pFac <- EnvFitRes[[BestMods[r, "mods"]]][[i]][[Ord]]$factors[["pvals"]]
      Values <- data.frame(r2 = c(RVec, RFac), p.value = c(pVec, pFac))
      Values <- signif(Values, digits = 3)
      Values <- cbind(var = rownames(Values), Values)
      rownames(Values) <- NULL
      Values <- cbind(mods = rep(BestMods[r, "mods"], nrow(Values)), Values)
      Values <- Values[order(Values[["var"]]),]
      ValuesRes <- rbind(ValuesRes, Values)
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
    ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] < 0.0001)] <-
      "<0.0001"
    ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] > 0.0001 &
                                   ValuesRes[["p.value"]] < 0.01)] <- "<0.01"
    ValuesRes[["p.value"]][which(ValuesRes[["p.value"]] > 0.01 &
                                   ValuesRes[["p.value"]] < 0.05)] <- "<0.05"
    return(ValuesRes)
  }), distances)
  VRList[[Ord]] <- Val
}
options(scipen = 0)

# generate the table (here one table per metric - e.g., richness, shannon, ... - including all normalisation methods)
if (length(list.dirs(file.path(DirSave, "Signif_Effects", "EnvFitBetaDivFunctionsBestModTab"))) == 0) {
  dir.create(file.path(DirSave, "Signif_Effects", "EnvFitBetaDivFunctionsBestModTab"))
} 

for(Ord in c("PCOA", "NMDS")){
  for(m in distances){
    Values <- VRList[[Ord]]
    ValIndex <- which(names(VRList[[Ord]]) == m)
    l <- vector("list", nrow(BestMods))
    l[[1]] <- kableExtra::kable_classic_2(
      kableExtra::kbl(Values[[ValIndex]][-1], 
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
    
    for (i in seq(nrow(BestMods))) {
      l[[i + 1]] <- kableExtra::pack_rows(l[[i]],
                                          BestMods[i, "modsId"],
                                          min(which(Values[[ValIndex]][[1]] ==  BestMods[i, "mods"])),
                                          max(which(Values[[ValIndex]][[1]] ==  BestMods[i, "mods"])),
                                          label_row_css = "background-color: #fff; color: #000; border-bottom: 1px solid; border-top: 1px solid; font-style: italic;")
    }
    l[[length(l)]]
    kableExtra::save_kable(l[[length(l)]], 
                           file.path(DirSave, "Signif_Effects", "EnvFitBetaDivFunctionsBestModTab", 
                                     paste0(paste("EnvFitBetaDivFunctionsBestModTab", m, Ord, sep = "_"), 
                                            ".HTML")))
  }
}

# generate the table including all metric for the original data log transformed (Table 5)
## retrieve the summary from the distances of the original Log transformed data in PCOA ordination (for both axes)
dataNorm <- "original.log"
ordination <- "PCOA"

ValuesPart <- VRList[[ordination]]
OrLog <- lapply(ValuesPart, function(x)
  x[which(x[[1]] == dataNorm), ])
OrLogDf <- do.call(rbind,OrLog)
OrLogDf$index <- gsub("\\.\\d+", "", rownames(OrLogDf))
rownames(OrLogDf) <- NULL

l <- vector("list", length(unique(OrLogDf[["index"]])))
l[[1]] <- kableExtra::kable_classic_2(
  kableExtra::kbl(
    OrLogDf[-c(1, length(OrLogDf))],
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
                       file.path(DirSave, "Signif_Effects", "EnvFitBetaDivFunctionsBestModTab", 
                                 paste0(paste("EnvFitBetaDivFunctionsBestModTab", dataNorm, ordination, sep = "_"), 
                                        ".HTML")))
