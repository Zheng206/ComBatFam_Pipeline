% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visual_prep.R
\name{visual_prep}
\alias{visual_prep}
\title{Batch Effect Visualization Datasets Preparation}
\usage{
visual_prep(
  type,
  features,
  batch,
  covariates,
  interaction = NULL,
  random = NULL,
  smooth = NULL,
  smooth_int_type = "linear",
  df,
  family = "comfam",
  ref.batch = NULL,
  predict = FALSE,
  object = NULL,
  reference = NULL,
  ...
)
}
\arguments{
\item{type}{A model function name that is used or to be used in the ComBatFamily Package (eg: "lmer", "lm").}

\item{features}{Features to be harmonized. \emph{n x p} data frame or matrix of observations where \emph{p} is the number of features and \emph{n} is the number of subjects.}

\item{batch}{Factor indicating batch (often equivalent to site or scanner).}

\item{covariates}{Name of covariates supplied to \code{model}.}

\item{interaction}{Expression of interaction terms supplied to \code{model} (eg: "age,diagnosis").}

\item{random}{Variable name of a random effect in linear mixed effect model.}

\item{smooth}{Variable name that requires a smooth function.}

\item{smooth_int_type}{Indicates the type of interaction in \code{gam} models. By default, smooth_int_type is set to be "linear", representing linear interaction terms.
"categorical-continuous", "factor-smooth" both represent categorical-continuous interactions ("factor-smooth" includes categorical variable as part of the smooth),
"tensor" represents interactions with different scales, and "smooth-smooth" represents interaction between smoothed variables.}

\item{df}{Dataset to be harmonized.}

\item{family}{combat family to use, comfam or covfam.}

\item{ref.batch}{reference batch.}

\item{predict}{A boolean variable indicating whether to run ComBat from scratch or apply existing model to new dataset (currently only work for original ComBat and ComBat-GAM).}

\item{object}{Existing ComBat model.}

\item{reference}{Dataset to be considered as the reference group.}

\item{...}{Additional arguments to \code{comfam} or \code{covfam} models.}
}
\value{
\code{visual_prep} returns a list containing the following components:
\item{summary_df}{Batch sample size summary}
\item{residual_add_df}{Residuals that might contain additive and multiplicative joint batch effects}
\item{residual_ml_df}{Residuals that might contain multiplicative batch effect}
\item{pr.feature}{PCA results}
\item{pca_df}{A dataframe contains features in the form of PCs}
\item{tsne_df}{A dataframe prepared for T-SNE plots}
\item{kr_test_df}{A dataframe contains Kenward-Roger(KR) test results}
\item{fk_test_df}{A dataframe contains Fligner-Killeen(FK) test results}
\item{mdmr.summary}{A dataframe contains MDMR results}
\item{anova_test_df}{A dataframe contains ANOVA test results}
\item{kw_test_df}{A dataframe contains Kruskal-Wallis test results}
\item{lv_test_df}{A dataframe contains Levene's test results}
\item{bl_test_df}{A dataframe contains Bartlett's test results}
\item{red}{A parameter to highlight significant p-values in result table}
\item{info}{A list contains input information like batch, covariates, df etc}
\item{eb_df}{A dataframe contains empirical Bayes estimates}
\item{com_family}{ComBat family to be considered: comfam, covfam}
\item{harmonized_df}{Harmonized dataset}
\item{combat.object}{saved ComBat model}
}
\description{
Prepare relevant datasets for batch/site effect visualization.
}
