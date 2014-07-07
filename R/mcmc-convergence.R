#' @title Computes convergence diagnistics
#' 
#' @description
#' Computes a very simple convergence diagnistic based on a t-test for equal means for different
#' windows of the chain. For each window, each parameter is compared with the samples from the last
#' half of the data and the t-test p-value is computed.
#' 
#' @details
#' Keep in mind that this function only compare means. It is not a complete diagnistic of convergence
#' since samples from two windows can be very different and still have the same mean. Generally, though,
#' the CoalHMM runs tends to have very different means if they are not converged.
#' 
#' @param  samples     Samples from a CoalHMM MCMC run.
#' @param  window_size Size of windows to use in the t-tests.
#' 
#' @return A data frame with p-values for the t-tests in each window.
#' 
#' @export
convergence_t_tests <- function(samples, window_size = nrow(samples)/10) {
  parameters <- parameter_samples(samples)
  
  no_windows <- floor(nrow(samples)/window_size/2)
  
  last_samples <- parameters[(nrow(parameters)/2):nrow(parameters),]
  p_values <- matrix(nrow=no_windows, ncol=ncol(parameters) + 2)
  colnames(p_values) <- c('window.begin', 'window.end', names(parameters))
  
  for (window in 1:no_windows) {
    begin <- window_size * (window - 1) + 1
    end <- begin + window_size - 1
    p_values[window,"window.begin"] <- begin
    p_values[window,"window.end"] <- end
    
    for (i in 1:ncol(parameters)) {
      p_values[window,i+2] <- t.test(parameters[begin:end,i], last_samples[,i])$p.value
    }
  }
  
  return(p_values)
}