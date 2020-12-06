#' unnormalized log Sigmoid Posterior
#'
#' Poserior with the following assumptions:
#' * logistic Regression modell, using sigmoid function to transform mapping from
#'  (-Inf, Inf) to [0, 1]
#' * Priors for Intercept and Slope Coefficients are weak and zero-mean normal distributed
#' with uncerrelated Covariance sigma
#'
#' @param position baring the coordinates of the posterior parameter
#' @param sigma double, diagonalelements for the covariace matrix indicating scale of the parameters priors
#' @inheritParams naive_nouturn_sampler
#'
#' @example
#' library(mlbench)
#' p <- mlbench.spirals(300,1.5,0.05)
#' plot(p)
#' design <- cbind(1, p$x)
#' target <- as.numeric(p$classes)
#' target[which(target == 2)] = 0
#' position <- rnorm(3)
#' @return unnormalized log posterior density value
sigmoid_posterior <- function(position, design = NULL, target = NULL, sigma = 200L) {
  if(is.data.frame(position)) position <- as.numeric(position)
  log_lik <- log(((1 + exp((-target) %*% (design %*% position))))^(-1))
  posterior <- (log_lik - ((2 * sigma ^ 2) ^(-1)) * t(position) %*% position)
  as.numeric(posterior)
}

#' @inheritParams sigmoid_posterior
#'
partial_deriv_sigmoid <- function(position, design = NULL, target = NULL, sigma = 200L) {
  if(is.data.frame(position)) position <- as.numeric(position)
  sigmoid_exp <- exp((target) %*% (design %*% position))
  sapply(seq_len(length(position)), function(x) {
    numerator <- (t(target) %*% design[,x])
    denominator <- 1 + sigmoid_exp
    penalizer <- position[x] / sigma
    as.numeric(numerator / denominator - penalizer)
  })
}
