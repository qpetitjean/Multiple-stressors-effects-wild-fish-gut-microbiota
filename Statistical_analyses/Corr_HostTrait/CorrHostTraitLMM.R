##############################################################################################################
#  Check Relationships with alteration of host’s traits                                                      #
#############################################################################################################
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/multiNorm.R") # import a function used to compute multiple normalization methods
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/DredgedLMM.R")  # import a function used to retrieve the summary, anova table and Rsquared of the most straightforward mixed effect model using rank selection (e.g., AICc)
source("D:/POSTDOC_INP_GOLFECH_2023/R scripts/GitRepository/R_Func/alphaDiv.R") # import a function compute alpha diversity indices from a matrix containing read counts per sample
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
rownames(labFish[["reads"]]) <- labFish$samples$Num_prlvt_Euth


#######################################
# compute multiple normalization      #
#######################################
normalizedDat <- multiNorm(
  dataList = labFish,
  norm = c(
    "original.log"
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

#############################################################
# Merge data with all host's traits and alpha div indices   #
#############################################################

FinalDat <- merge(
  x = FullDat,
  y = alphaDivindices[["original.log"]],
  by.x = "Ind",
  by.y = "Samples",
  all = TRUE
)

# split the dataset between C and NC contamination treatment
aDivC <- FinalDat[which(FinalDat$Contam == "C"), ]
aDivNC <- FinalDat[which(FinalDat$Contam == "NC"), ]

#############################################################################################################
# Testing correlation between host traits and Alpha div among C and NC treatments                           #
#############################################################################################################
# for more information about host's traits measured in this experiment
# see Petitjean et al., 2021, Intraspecific variability of responses to combined metal contamination and immune challenge among wild fish populations, Environmental Pollution, Volume 272, 2021, 116042, ISSN 0269-7491, https://doi.org/10.1016/j.envpol.2020.116042.

#########################
# for fish activity - pca axis explained mostly by the time spent swimming 
Traits <-
  c("General_activity_Axis1_after_treat",
    "Voracity_Axis2_after_treat",
    "Sociability_Axis3_after_treat",
    "Movement_second_after_treat",
    "ANND_cm_after_treat",
    "Latency_time_second_after_treat",
    "Nbr_Eaten_Pellets_after_treat",
    "ratio_H_L",
    "nb_lymphocytes",
    "nb_neutrophiles",
    "AE_mj.mg",
    "Lipides_mj.mg",
    "Prot_mj.mg",
    "Carbo_mj.mg",
    "DailyMLoss",
    "H2O2_mM", 
    "mMHCLO")

CorrRes <- list()
for(HostT in Traits){
for(i in indices){
  for(treat in c("NC", "C")){
    dataTemp <- get(paste0("aDiv", treat))
    dataTemp <- dataTemp[!is.na(dataTemp[[i]]), ]
    dataTemp <- dataTemp[!is.na(dataTemp[[HostT]]), ]
    dataTemp[["Bac"]] <- as.character(dataTemp[["Bac"]])
    
    ModTemp <- lme4::lmer(
      as.formula(paste0(
        paste0(HostT, " ~ "), i, paste0("+", "(1 |Session / Bac) + (1 | Pop)")
      )),
      data = dataTemp,
      na.action = "na.fail",
      REML = T
    )
    
    ##### create output folder to save figs
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota"))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota"))
    }
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota", HostT))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota", HostT))
    }
    if (length(list.dirs(file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota", HostT, i))) == 0) {
      dir.create(file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota", HostT, i))
    }
    toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota", HostT, i, paste0("Mod_Perf_", treat, ".tif" ))
    if(file.exists(toSave)){
      unlink(toSave)
    } 
    grDevices::tiff(filename = toSave,
                    width = 1200,
                    height = 1200,
                    res = 100,
                    compression = "lzw",
                    pointsize = 8)
    print(performance::check_model(ModTemp))
    dev.off()
    
    
    AOV <- car::Anova(ModTemp, type = "3")
    Sum <- summary(ModTemp)
    # retrieve the results in a list
    CorrRes[[HostT]][[i]][[treat]] <- list(
      SampleSize = Sum$devcomp$dims[["n"]],
      ModL = ModTemp,
      ModSum = Sum,
      ModAnov = AOV,
      Rsquared = MuMIn::r.squaredGLMM(ModTemp)
    )
    
    toSave <- file.path(DirSave, "Signif_Effects", "Signif_Plots", "CorrHostMicrobiota", HostT, i, paste0(paste("Corr", treat, sep = "_"), ".tif" ))
    if(file.exists(toSave)){
      unlink(toSave)
    } 
    grDevices::tiff(filename = toSave,
                    width = 1200,
                    height = 800,
                    res = 250,
                    compression = "lzw",
                    pointsize = 10)
    plot(dataTemp[[HostT]]~ dataTemp[[i]],
         xlab = i,
         ylab = HostT,
         pch = 19, 
         las = 1)
    mtext(paste("treatment:", treat), side = 3, line = 2)
    p <- paste("p-value=", signif(AOV$`Pr(>Chisq)`[2], 2))
    est <- paste("Estimate=", signif(Sum$coefficients[grep(i, rownames(Sum$coefficients)), "Estimate"], 2))
    mtext(paste(p, est, sep = " "), side = 3, )
    if(AOV$`Pr(>Chisq)`[2] <0.05){
    abline(coef = c( 
      Sum$coefficients["(Intercept)", "Estimate"],
      Sum$coefficients[grep(i, rownames(Sum$coefficients)), "Estimate"]), 
      col = "firebrick")
    }
    dev.off()
  }
}
}
