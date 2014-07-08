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
convergence_t_tests <- function(samples, window_size = nrow(samples)/20) {
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

#' @title Determines if the diagnistics means that we reject convergence.
#' 
#' @description
#' Since the t-test diagnostics give a p-value per parameter we use this function to determine
#' if the parameters combined leads us to reject convergence. This is simply done by correcting
#' for multiple tests and then using the smallest p-value per row.
#' 
#' @param convergence_p_values  P-values computed as in \code{convergence_t_tests}.
#' @return A data frame with a row per window specifying if convergence was rejected.
#' 
#' @export
rejected_convergence <- function(convergence_p_values) {
  
  # Bonferroni correction for multiple tests
  no_parameters <- ncol(convergence_p_values) - 2
  adjusted_p_values <- convergence_p_values[,-c(1,2)] * no_parameters
  
  data.frame(convergence_p_values[,c(1,2)],
             rejected = apply(adjusted_p_values, 1, function(x) any(x < 0.05)))
}

#' @title Computes the probability of convergence at each windows
#' 
#' @description
#' Computes the probability, at each window, that the chain was converged at that point, assuming
#' that the \emph{has} converged at some window.
#' 
#' @details
#' Computes the probability that convergence happened at a given point based on the pattern of 
#' rejected convergence tests. It assumes that the rejection is picked with a probability that is
#' beta(a,b) distributed before the change point and rejected with 5% probability after the change
#' point.
#' 
#' @export
convergence_probabilities <- function(rejections, a=1, b=1) {
  no_points <- length(rejections)
  probabilities <- numeric(no_points)
  probabilities[1] <- prod(0.05**rejections)
  for (point in 2:no_points) {
    before <- rejections[1:(point-1)]
    after <- rejections[point:no_points]
    probabilities[point] <- prod(beta(a + before, b + (1 - before))) * prod(0.05**after)
  }
  not_converged <- prod(beta(a + rejections, b + (1 - rejections)))
  
  total_probability <- sum(probabilities) + not_converged
  
  result <- list(change_probabilities = probabilities / total_probability,
                 not_converged = not_converged / total_probability)
  
  return(result)
}

#' @title Computes convergence diagnostics
#' 
#' @description
#' Computes convergence diagnistics based on the \code{convergence_t_tests} diagnostics and finds the
#' most likely convergence point using \code{convergence_probabilities}.
#' 
#' @param samples        Samples from a CoalHMM MCMC run.
#' @param window_size    Window size used for the convergence diagnostics.
#' @param a              Meta-parameter for the prior probability of rejection before convergence.
#' @param b              Meta-parameter for the prior probability of rejection before convergence.
#' 
#' @export
convergence_diagnostics <- function(samples, window_size = nrow(samples)/20, a = 1, b = 1) {
  p_vals <- convergence_t_tests(samples, window_size)
  converged <- rejected_convergence(p_vals)
  change_probs <- convergence_probabilities(converged$rejected)
  converged$change_point_probabilities <- change_probs$change_probabilities
  change_point <- converged$window.end[which.max(converged$change_point_probabilities)]
    
  result <- list(convergence = converged, convergence_point = change_point,
                 not_converged = change_probs$not_converged)
  class(result) <- "convergence_diag"
  
  return(result)
}

plot.convergence_diag <- function(diag, ...) {
  plot(diag$convergence$window.end, diag$convergence$change_point_probabilities,
       main='Convergence point probabilities',
       xlab='Convergence point', ylab='Convergence probability',
       type='o', pch=20,
       ...)
  abline(v=diag$convergence_point, col='red', lty='dashed')
}

#' @title Remove the estimated burn-in samples
#' 
#' @description
#' Estimate the convergence point of the chain and remove the burn-in samples.
#' 
#' @export
strip_burnin <- function(samples, ...) {
  diag <- convergence_diagnostics(samples, ...)
  return(samples[diag$convergence_point:nrow(samples),])
}
