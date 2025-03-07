#' Retrieve Model Summary, ANOVA, and R-squared via Model Dredging
#'
#' @description
#' \code{DredgedLMM} fits a full mixed effects model using \code{lmer} from the \pkg{lme4} package and
#' then performs model selection using the \code{dredge} function from the \pkg{MuMIn} package. Depending on the
#' specified method, it selects the final model among those with a delta value below 2 either by the highest weight ("weight")
#' or by the lowest ranking criterion ("parsimony"). The function refits the model using the selected predictors and returns a list
#' containing the sample sizes, the final model, its summary, a Type III ANOVA table (via \pkg{car}), and the marginal and conditional
#' R-squared values.
#'
#' @param Data A data frame containing the dataset to be modeled.
#' @param expVar A character string specifying the response variable.
#' @param respExpr A character string specifying the fixed-effects portion of the model (including main effects and interaction terms).
#' @param random A character string specifying the random effects structure (e.g., "(1 | Session/Bac) + (1 | Pop)").
#' @param rank A character string specifying the ranking criterion for model selection (default is "AICc").
#' @param method A character string indicating the method for final model selection when several candidate models have delta values below 2.
#'        Options are "weight" (select the model with the highest weight) or "parsimony" (select the model with the lowest rank value).
#'
#' @return A list with the following components:
#' \describe{
#'   \item{SampleSize}{A data frame reporting sample sizes for each predictor variable included in the final model.}
#'   \item{ModL}{The final refined mixed effects model (an object of class \code{lmerMod}).}
#'   \item{ModSum}{The summary of the final model.}
#'   \item{ModAnov}{The Type III ANOVA table for the final model.}
#'   \item{Rsquared}{A numeric vector containing the marginal and conditional R-squared values as computed by \code{MuMIn::r.squaredGLMM}.}
#' }
#'
#' @details
#' \code{DredgedLMM} first filters the input data to exclude missing values in the response variable.
#' A full model is then fitted using the formula constructed from the response variable, the specified fixed effects, and random effects.
#' Model selection is performed with \code{MuMIn::dredge}, and if multiple models have a delta ranking below 2,
#' the final model is chosen based on the selected method ("weight" or "parsimony"). The chosen predictors are then used
#' to refit the model using REML = TRUE, and the function returns the refined model along with its summary, ANOVA table,
#' and R-squared values.
#'
#' **Dependencies:**
#' This function depends on several packages. It internally calls functions from:
#' \itemize{
#'   \item \strong{lme4} (for linear mixed model - LMM - model fit),
#'   \item \strong{MuMIn} (for model selection and r2 computation),
#'   \item \strong{car} (for p-value computation from LMM)
#' }
#'
#' For reproducibility, please ensure that the required packages are installed and loaded.
#' Although not ideal to load packages within a function, the following libraries are explicitly 
#' loaded here to ensure that all dependent functions are available:
#'
#' @importFrom lme4 lmer
#' @importFrom MuMIn dredge r.squaredGLMM
#' @importFrom car Anova 
#' 
#' @author Quentin PETITJEAN [quentin.petitjean@inrae.fr]
#'
#' @date 15/06/2023
#' 
#' @export
 
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
