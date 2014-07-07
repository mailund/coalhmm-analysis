
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
  samples[,!(names(samples) %in% c("sample", "log.prior","log.likelihood","log.posterior"))]
}

#' @title Creates an object for manipulating CoalHMM MCMC chains
#' 
#' @param samples   Samples from a CoalHMM MCMC run.
#' @return The samples wrapped as an object of class "coalhmm_chain".
#' @export
coalhmm_chain <- function(samples) {
  samples$sample <- 1:nrow(samples)
  class(samples) <- c('coalhmm_chain', 'data.frame')
  return(samples)
}

plot.coalhmm_chain <- function(chain, ...) {
  convergence_diag <- convergence_diagnostics(chain, ...)
  convergence_point <- convergence_diag$convergence_point
  parameters <- parameter_samples(chain)
  melted <- melt(parameters, id='sample')
  ggplot(melted, aes(x=sample, y=value)) + geom_point() + 
    facet_grid(variable~., scales = "free_y") + geom_smooth(method="loess") +
    geom_vline(xintercept=convergence_point, color='red') + theme_bw()
}

summary.coalhmm_chain <- function(chain, ...) {
  convergence_diag <- convergence_diagnostics(chain, ...)
  convergence_point <- convergence_diag$convergence_point
  parameters <- parameter_samples(chain[convergence_point:nrow(chain), ])

  estimates <- apply(parameters, 2, quantile, probs=c(0.025, 0.25, 0.5, 0.75, 0.975))
  
  result <- list(parameters = parameters, estimates = estimates,
                 convergence_point = convergence_point)
  class(result) <- "coalhmm_chain_summary"
  return(result)
}

print.coalhmm_chain_summary <- function(chain_summary) {
  cat("Estimated convergence point is", chain_summary$convergence_point)
  cat(" leaving", nrow(chain_summary$parameters), "samples\n\n")
  
  parameter_names <- colnames(chain_summary$estimates)
  cat("Estimates (quantiles):\n")
  print(chain_summary$estimates)
  cat("\n")
}

plot.coalhmm_chain_summary <- function(chain_summary) {
  ggplot(data=melt(chain_summary$parameters)) + geom_histogram(aes(x=value)) +
    facet_grid(. ~ variable, scale="free") + theme_bw()
}