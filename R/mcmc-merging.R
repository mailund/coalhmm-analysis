#' @title Compares the distributions for each parameter for two chains to determine if they are equal
#' 
#' @description
#' Compares the samples of two chains to test that for each parameter the samples are from the same distribution.
#' 
#' @details
#' For each parameter the test does a Kolmogorov-Smirnov test (\code{ks.test}). This is a very simple approach to testing
#' equality as it doesn't take into account the dependencies between parameters nor that the samples are not independent,
#' but it is a simple and fast test. It should really be supplemented with visual inspection and other diagnostics.
#' 
#' @param chain1         Samples from an MCMC run.
#' @param chain2         Samples from an MCMC run.
#' @param strip_chains   Should the chains be stripped of burn-in samples before comparison?
#' 
#' @export
equal_distribution <- function(chain1, chain2, strip_chains = TRUE, ...) {
  if (strip_chains) {
    parameters1 <- parameter_samples(chain1)
    parameters2 <- parameter_samples(chain2)
    
  } else {
    parameters1 <- parameter_samples(strip_burnin(chain1, ...))
    parameters2 <- parameter_samples(strip_burnin(chain2, ...))
  }
  
  # FIXME: check that they actually have the same parameters
  
  param_names <- names(parameters1)
  p_values <- numeric(length(param_names))
  for (p in 1:length(param_names)) {
    p_values[p] <- ks.test(parameters1[,p], parameters2[,p], ...)$p.value
  }
  names(p_values) <- param_names
  result <- list(p_values = p_values)
  class(result) <- "equal_distribution_test"
  return(result)
}

print.equal_distribution_test <- function(test) {
  print(test$p_values)
}

summary.equal_distribution_test <- function(test) {
  class(test) <- "equal_distribution_test_summary"
  return(test)
}

print.equal_distribution_test_summary <- function(test) {
  cat("Individual tests, with non-adjusted p-values:\n\n")
  print(test$p_values)
  cat("\n")
}

#' @title Tests if a list of chains can be merged
#' 
#' @description
#' Test, for all pairs of chains, if they are samples from the same distribution (using \code{equal_distribution}).
#' 
#' @details
#' The return value contains a p-value computed from another \code{ks.test} testing if the test p-values are uniform.
#' The idea here is that if each parameter is actually equally distributed then the p-values should be uniform.
#' 
#' @param chains  A \code{list} containing \code{coalhmm_chain}s.
mergeable_chains <- function(chain_list, strip_chains=TRUE, ...) {
  N <- length(chain_list)
  parameters <- colnames(parameter_samples(chain_list[[1]]))
  p_values <- matrix(nrow = N * (N - 1) / 2, ncol = length(parameters) + 2)
  row <- 1
  for (chain1 in 1:(N-1)) {
    for (chain2 in (chain1+1):N) {
      p_values[row,] <- c(chain1, chain2, 
                          equal_distribution(chain_list[[chain1]], chain_list[[chain2]], strip_chains, ...)$p_values)
      row <- row + 1
    }
  }
  colnames(p_values) <- c("chain_1", "chain_2", parameters)
  result <- list(p_values = p_values)
  class(result) <- "mergeable_chains_test"
  return(result)
}

print.mergeable_chains_test <- function(test) {
  print(test$p_values)
}

merge_chains <- function(chain_list, strip_chains = TRUE) {
  if (strip_chains) {
    return(do.call(rbind, lapply(chain_list, strip_burnin)))
  } else {
    return(do.call(rbind, chain_list))
  }
}