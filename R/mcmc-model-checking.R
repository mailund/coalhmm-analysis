
#' @title Computes the model log likelihood from MCMC posterior samples
#' 
#' @description
#' Computes the so-called model likelihood, that is the probability of the data given the model,
#' \eqn{P(D | M) = \int P(D | \theta, M) P(\theta | M) d\theta}, based on samples from an MCMC run.
#' The ratio between two model likelihoods from different models is the Bayes factor and can be
#' used for model selection.
#'
#' @details
#' Computes the so-called model likelihood, that is the probability of the data given the model,
#' \eqn{P(D | M) = \int P(D | \theta, M) P(\theta | M) d\theta}.
#' 
#' Since the samples from the MCMC are from the posterior and not the prior distribution the likelihood
#' cannot be computed simply as the mean of the likelihood -- and typically that would be a very
#' inefficient approach to computing the likelihood -- so instead the harmonic mean of the likelihood
#' sampled from the posterior is used. This follows from 
#' \deqn{\int P(\theta | D, M) / P(D | \theta, M) d\theta 
#'  = \int 1/P(D | M) [P(D | \theta, M)P(\theta | M)]/P(D | \theta, M) d\theta
#'  = 1/P(D | M)
#' } 
#' so \eqn{1/P(D | M) = \int 1/P(D | \theta, M) P(\theta | D, M) d\theta}
#' which is approximated from samples from the posterior \eqn{P(\theta | D, M)}.
#' 
#' @param samples   A data frame with the samples from an MCMC run.
#' @param burn.in   The number of samples to discard as burn-in samples.
#' 
#' @return The model log-likelihood \eqn{\log P(D | M)}.
#' 
#' @export
model.likelihood <- function(samples, burn.in = 0) {
  
  # Discard burn-in samples.
  log.likelihoods <- samples$log.likelihood[(1+burn.in):nrow(samples)]
  N <- length(log.likelihoods)
  
  scale <- -mean(log.likelihoods)
  scaled.sum <- sum(1/exp(log.likelihoods + scale))
  log.harmonic.mean <- scale - log(N) + log(scaled.sum)
  
  # P(D | M) = 1/"harmonic mean" so log(P(D | M)) = log(1) - log("harmonic mean")
  -log.harmonic.mean
}