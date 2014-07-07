
#' @title Extract the parameters from the samples data frame
#' 
#' @description
#' Pulls out the columns that corresponds sampled parameters rather than the probabilities,
#' that is it returns a data frame without the "log.prior", "log.likelihood" and "log.posterior"
#' columns.
#' 
#' @details
#' For packages that implements convergence diagnostics, such as coda, we need the samples with
#' only the parameters and not the prior and posterior probabilities or the likelihood.
#' 
#' This function simply returns that from a data frame of samples as output by the CoalHMM MCMC
#' samplers.
#' 
#' @param samples  Data frame with samples from a CoalHMM MCMC.
#' @return The samples containing only the parameter columns.
#' @export
parameter_samples <- function(samples) {
  samples[,!(names(samples) %in% c("log.prior","log.likelihood","log.posterior"))]
}

