% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Harmonization.R
\name{combat_harm}
\alias{combat_harm}
\title{ComBatFamily Harmonization}
\usage{
combat_harm(
  result = NULL,
  features = NULL,
  batch = NULL,
  covariates = NULL,
  df = NULL,
  type = "lm",
  random = NULL,
  smooth = NULL,
  interaction = NULL,
  smooth_int_type = NULL,
  family = "comfam",
  ref.batch = NULL,
  predict = FALSE,
  object = NULL,
  reference = NULL,
  out_ref_include = TRUE,
  ...
)
}
\arguments{
\item{result}{A list derived from \code{visual_prep()} that contains datasets for shiny visualization. Can be skipped if \code{features}, \code{batch}, \code{covariates} and \code{df} are provided.}

\item{features}{Features to be harmonized. \emph{n x p} data frame or matrix of observations where \emph{p} is the number of features and \emph{n} is the number of subjects. Can be skipped if \code{result} is provided.}

\item{batch}{Factor indicating batch (often equivalent to site or scanner). Can be skipped if \code{result} is provided.}

\item{covariates}{Name of covariates supplied to \code{model}. Can be skipped if \code{result} is provided.}

\item{df}{Dataset to be harmonized. Can be skipped if \code{result} is provided.}

\item{type}{A model function name that is used or to be used in the ComBatFamily Package (eg: "lmer", "lm").}

\item{random}{Variable name of a random effect in linear mixed effect model.}

\item{smooth}{Variable name that requires a smooth function.}

\item{interaction}{Expression of interaction terms supplied to \code{model} (eg: "age,diagnosis").}

\item{smooth_int_type}{A vector that indicates the types of interaction in \code{gam} models. By default, smooth_int_type is set to be NULL, "linear" represents linear interaction terms.
"categorical-continuous", "factor-smooth" both represent categorical-continuous interactions ("factor-smooth" includes categorical variable as part of the smooth),
"tensor" represents interactions with different scales, and "smooth-smooth" represents interaction between smoothed variables.}

\item{family}{combat family to use, comfam or covfam.}

\item{ref.batch}{reference batch.}

\item{predict}{A boolean variable indicating whether to run ComBat from scratch or apply existing model to new dataset (currently only work for original ComBat and ComBat-GAM).}

\item{object}{Existing ComBat model.}

\item{reference}{Dataset to be considered as the reference group.}

\item{out_ref_include}{A boolean variable indicating whether the reference data should be included in the harmonized data output.}

\item{...}{Additional arguments to \code{comfam} or \code{covfam} models.}
}
\value{
\code{combat_harm} returns a list containing the following components:
\item{eb_df}{A dataframe contains empirical Bayes estimates}
\item{com_family}{ComBat family to be considered: comfam, covfam}
\item{harmonized_df}{Harmonized dataset}
\item{combat.object}{saved ComBat model}
}
\description{
Conduct harmonization using the ComBatFamily package, which currently includes 4 types of harmonization techniques: 1) Original ComBat, 2) Longitudinal ComBat, 3) ComBat-GAM, 4) CovBat.
}
