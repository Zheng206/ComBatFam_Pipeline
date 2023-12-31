% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/visual_prep.R
\name{form_gen}
\alias{form_gen}
\title{ComBatFam Formula generations}
\usage{
form_gen(x, c, i, random, smooth, int_type)
}
\arguments{
\item{x}{A model function name that is used or to be used in the ComBatFamily Package (eg: "lmer", "lm").}

\item{c}{Data frame or matrix of covariates supplied to \code{model}}

\item{i}{Expression of interaction terms supplied to \code{model}, using comma to separate terms. (eg: "age,diagnosis").}

\item{random}{Variable name of a random effect in linear mixed effect model.}

\item{smooth}{Variable name that requires a smooth function.}

\item{int_type}{Indicates the type of interaction in \code{gam} models. By default, smooth_int_type is set to be "linear", representing linear interaction terms. "factor-smooth" represents categorical-continuous interactions and "tensor" represents interactions with different scales.}
}
\value{
A string of formula
}
\description{
Generate appropriate formula for ComBatFam models
}
