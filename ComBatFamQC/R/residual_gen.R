require(stats)
require(lme4)
require(parallel)
require(mgcv)
require(dplyr)
require(tidyverse)
require(stringr)

#' Post Harmonization Residual Generation
#'
#' Extract residuals after harmonization.
#'
#' @param type A model function name that is to be used (eg: "lmer", "lm").
#' @param features Features/rois to extract residuals from. \emph{n x p} data frame or matrix of observations where \emph{p} is the number of features and \emph{n} is the number of subjects.
#' @param batch Factor indicating batch (often equivalent to site or scanner).
#' @param covariates Name of covariates supplied to `model`.
#' @param interaction Expression of interaction terms supplied to `model` (eg: "age*diagnosis").
#' @param random Variable name of a random effect in linear mixed effect model.
#' @param smooth Variable name that requires a smooth function.
#' @param smooth_int_type Indicates the type of interaction in `gam` models. By default, smooth_int_type is set to be "linear", representing linear interaction terms. "factor-smooth" represents categorical-continuous interactions and "tensor" represents interactions with different scales.
#' @param df Harmonized dataset to extract residuals from.
#' @param rm variables to remove effects from (apart from batch).
#' @param model A boolean variable indicating whether an existing model is to be used.
#' @param model_path path to the existing model.
#'
#' @return `residual_gen` returns a list containing the following components:
#' \item{model}{a list of regression models for all rois}
#' \item{residual}{Residual dataframe}
#' 
#' 
#' @import parallel
#' @import tidyverse
#' @import dplyr
#' @import stringr
#' @importFrom broom tidy
#' @importFrom mgcv gam 
#' @importFrom lme4 lmer
#' @importFrom stats lm median model.matrix prcomp predict qnorm update var anova as.formula coef resid na.omit complete.cases
#' 
#' @export
#' 
#' 

residual_gen <- function(type, features, batch, covariates, interaction = NULL, random = NULL, smooth = NULL, smooth_int_type = "linear", df, rm = NULL, model = FALSE, model_path = NULL){
  obs_n = nrow(df)
  df = df[complete.cases(df[c(features, batch, covariates, random)]),]
  obs_new = nrow(df)
  print(paste0(obs_n - obs_new, " observations are dropped due to missing values.")) 
  df[[batch]] = as.factor(df[[batch]])
  char_var = covariates[sapply(df[covariates], function(col) is.character(col) || is.factor(col))]
  enco_var = covariates[sapply(df[covariates], function(col) length(unique(col)) == 2 && all(unique(col) %in% c(0,1)))]
  df[char_var] =  lapply(df[char_var], as.factor)
  df[enco_var] =  lapply(df[enco_var], as.factor)
  cov_shiny = covariates
  char_var = c(char_var, enco_var)
  if(!is.null(random)){
    for (r in random){
      df[[r]] = as.factor(df[[r]])
    }
  }
  
  if(!is.null(interaction)){
    if(type == "lmer"){
      #df[random] = sapply(random, function(x) as.factor(df[[x]]), USE.NAMES = FALSE)
      interaction = gsub(",", ":", interaction)
    }else if(type == "gam"){
      covariates = setdiff(covariates, smooth)
      if(smooth_int_type == "linear"){
        interaction = gsub(",", ":", interaction)
      }else if(smooth_int_type == "categorical-continuous"){
        element = lapply(interaction, function(x) str_split(x,",")[[1]])
        smooth_element = lapply(element, function(x) x[which(x %in% smooth)])
        categorical_element = lapply(1:length(element), function(i) setdiff(element[[i]], smooth_element[[i]]))
        interaction = lapply(1:length(element), function(i) paste0("s(", smooth_element[[i]], ", by = ", categorical_element[[i]], ")")) %>% unlist()
        smooth = setdiff(smooth, unique(unlist(smooth_element)))
      }else if(smooth_int_type == "factor-smooth"){
        interaction = paste("s(", interaction, ", bs = 'fs')")
        smooth_index = sapply(smooth, function(x){
          if(sum(grepl(x, interaction)) > 0){return(1)}else(return(0))
        }, USE.NAMES = FALSE)
        cov_index = sapply(covariates, function(x){
          if(sum(grepl(x, interaction)) > 0){return(1)}else(return(0))
        }, USE.NAMES = FALSE)
        smooth = smooth[which(smooth_index == 0)]
        covariates = covariates[which(cov_index == 0)]
      }else if(smooth_int_type == "tensor"){
        interaction = paste("ti(", interaction, ")")
      }
    }else if(type == "lm"){interaction = gsub(",", ":", interaction)}
  }
  
  used_col = c(features)
  other_col = setdiff(colnames(df), used_col)
  other_info = df[other_col]
  rm = c(batch, rm)
  if(!model){
    models = mclapply(features, function(y){
      model = model_gen(y = y, type = type, batch = batch, covariates = covariates, interaction = interaction, random = random, smooth = smooth, df = df)
      return(model)
    }, mc.cores = detectCores())
  }else{
    models = readRDS(model_path)
  }
  
  if(type!="lmer"){
    residuals = mclapply(1:length(features), function(i){
      #sub_df = df %>% dplyr::select(features[i], batch, covariates)
      #for (x in rm){
      #  sub_df[x] = rep(0, nrow(sub_df))
      #}
      model_coef = coef(models[[i]])
      rm_names = c()
      for (x in rm){
        sub_name = names(model_coef)[which(grepl(x, names(model_coef)))]
        rm_names = c(rm_names, sub_name)
      }
      use_coef = model_coef[!names(model_coef) %in% rm_names]
      predict_y = model.matrix(models[[i]])[, which(grepl(paste0(names(use_coef), collapse = "|"), colnames(model.matrix(models[[i]]))))] %*% t(t(unname(use_coef)))
      #predict_y = predict(models[[i]], newdata = sub_df, type = "response")
      residual_y = df[[features[i]]] - predict_y
      residual_y = data.frame(residual_y)
    }, mc.cores = detectCores()) %>% bind_cols()
  }else{
    df[[random]] = as.factor(df[[random]])
    residuals = mclapply(1:length(features), function(i){
      model_coef = coef(models[[i]])[[1]]
      int = model_coef[1]
      random_effect = cbind(rownames(int), int)
      colnames(random_effect) = c(random, "intercept")
      random_effect[[random]] = as.factor(random_effect[[random]])
      model_coef = distinct(model_coef[-1])
      rm_names = c()
      for (x in rm){
        sub_name = names(model_coef)[which(grepl(x, names(model_coef)))]
        rm_names = c(rm_names, sub_name)
      }
      use_coef = model_coef[!names(model_coef) %in% rm_names]
      if(length(use_coef) == 0){predict_y = df[random] %>% left_join(random_effect, by = random) %>% pull(intercept)}else{
        predict_y = model.matrix(models[[i]])[, which(grepl(paste0(names(use_coef), collapse = "|"), colnames(model.matrix(models[[i]]))))] %*% t(use_coef)
      }
      #predict_y = predict(models[[i]], newdata = sub_df, type = "response")
      residual_y = df[[features[i]]] - predict_y
      residual_y = data.frame(residual_y)
    }, mc.cores = detectCores()) %>% bind_cols()
  }
  colnames(residuals) = features
  residuals = cbind(other_info, residuals)
  residuals = residuals[colnames(df)]
  result = list("model" = models, "residual"= residuals)
  return(result)
}

utils::globalVariables(c("features", "batch", "covariates", "intercept", "random"))
