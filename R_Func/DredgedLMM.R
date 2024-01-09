#############################################################################################################################################################################################
# This is function aims to retrieve the summary, anova table and Rsquared of the most straightforward mixed effect model using rank selection (e.g., AICc)                                  #
# Title: 
# author:  "Quentin PETITJEAN[q.petitjean1@gmail.com]                                                                                                                                       #
# date: "15/06/2023"                                                                                                                                                                        #
#############################################################################################################################################################################################

DredgedLMM <- function(Data = NULL, # a dataframe containing the data to model
                       expVar = NULL, # the explainatory variable
                       respExpr = NULL, # the full set of response variables, including interactions terms
                       random = NULL, # the random effect 
                       rank = "AICc", # see ?MuMIn::dredge 
                       method = c("weight", "parsimony")) { # the method to select the final model after dredging. If several model have a rank value below 2, the parsimony method return the model with the lowest rank value while the weight method returns the model with the highest weight 
  
  toMod <- Data[!is.na(Data[[expVar]]), ]
  
  # Build a full model and refined it using AICc
  Mod1 <-
    lme4::lmer(
      as.formula(paste0(
        paste0(expVar, " ~ "), respExpr, paste0("+", random)
      )),
      data = toMod,
      na.action = "na.fail",
      REML = F
    )
  
Dredged <- tryCatch({
  suppressWarnings(
    suppressMessages(
      # Using do.call to execute test_function
      do.call(MuMIn::dredge, c(Mod1,  rank = rank))
    )
  )
}, warning = function(w) {
  # Here you can handle warnings, or simply ignore them.
  # Returning NULL or any other value replaces the actual output of the do.call in case of a warning.
  return(NULL)
}, error = function(e) {
  # Here you can handle errors, or simply ignore them.
  # Returning NULL or any other value replaces the actual output of the do.call in case of an error.
  return(NULL)
})

  # in case two or more models have delta AICc below 2, select only one according to the method argument
  if (length(which(Dredged[["delta"]] < 2)) == 1) {
    df <- as.data.frame(Dredged[1])
  } else{
    if (method == "weight") {
      df <-
        as.data.frame(Dredged[which(Dredged[["delta"]] < 2)][which.max(Dredged[["weight"]])])
    } else if (method == "parsimony") {
      df <-
        as.data.frame(Dredged[which(Dredged[["delta"]] < 2)][which.min(Dredged[[rank]])])
    }
  }
  df <- df[, colSums(is.na(df)) == 0]
  toKeep <- names(df[,-which(names(df) %in% c("(Intercept)", "df", "logLik", rank, "delta", "weight"))])
  
  # retrieve sample size for each variable remaining in the final Model
  sampleSize <- c()
  for(i in seq_along(toKeep)){
    if(toKeep[[i]] %in% names(toMod)){
      keepVar <- toKeep[[i]]
    }else{
      keepVar <-  gsub(":", ".", toKeep[[i]])
    }
    tab <- as.data.frame(table(toMod[[keepVar]]))
    tab[["treat"]] <- rep(keepVar, nrow(tab))
      sampleSize <- rbind(sampleSize, tab) 
  }

  # build the refined model
  ModFinal <-
    lme4::lmer(
      as.formula(paste0(
        paste0(expVar, " ~ "),
        do.call("paste", c(as.list(toKeep), sep = " + ")),
        paste0("+", random)
      )),
      data = toMod,
      na.action = "na.fail",
      REML = T
    )
  
  # retrieve the results in a list
  Res <- list(
    SampleSize = sampleSize,
    ModL = ModFinal,
    ModSum = summary(ModFinal),
    ModAnov = car::Anova(ModFinal, type = "3"),
    Rsquared = MuMIn::r.squaredGLMM(ModFinal)
  )
  
  return(Res)
}
