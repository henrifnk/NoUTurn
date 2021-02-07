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
#'
sigmoid_posterior <- function(position, design = NULL, target = NULL, sigma = 200L) {
  if(is.data.frame(position)) position <- as.numeric(position)
  log_lik <- -sum(log(((1 + exp((-target) * (design %*% position))))))
  posterior <- (log_lik - ((2 * sigma) ^(-1)) * t(position) %*% position)
  as.numeric(posterior)
}

#' @inheritParams sigmoid_posterior
#'
partial_deriv_sigmoid <- function(position, design = NULL, target = NULL, sigma = 100L) {
  if(is.data.frame(position)) position <- as.numeric(position)
  sigmoid_exp <- 1 + exp((target) * (design %*% position))
  sapply(seq_len(length(position)), function(x) {
    numerator <- as.numeric(target * design[,x])
    sum <- -sum(numerator / sigmoid_exp)
    penalizer <- position[x] / sigma
    as.numeric(sum - penalizer)
  })
}

#' @see https://cran.r-project.org/web/packages/hmclearn/vignettes/logistic_regression_hmclearn.html
bernoulli_pen <- function(position, sigma = 1e4) {
  if(is.data.frame(position)) position <- as.numeric(position)
  comp1 <- as.numeric(t(position) %*% t(design) %*% (target - 1))
  comp2 <- sum(log(1 + exp(-(design %*% position))))
  penal <- as.numeric(((2 * sigma) ^(-1)) * t(position) %*% position)
  comp1 - comp2  - penal
}

partial_deriv_benoulli <- function(position, sigma = 1e4){
  sigmoid_exp <- exp(-(design %*% position))
  comp1 <- (target - 1) + (sigmoid_exp / (1 + sigmoid_exp))
  comp2 <- t(design) %*% comp1
  comp2 - (position / sigma)
}

gamma_distr <- function(position, design, target, sigma = 100L){
  if(is.data.frame(design)) design <- as.matrix(design)
  dispersion <- position[1]
  beta <- position[-1]
  lin_predict <- design %*% beta
  comp1 <- (-dispersion) * (sum(target * exp(-lin_predict)) + sum(lin_predict))
  penal <- as.numeric(((2 * sigma) ^(-1)) * t(position) %*% position)
  sum(comp1) - penal
}

partial_deriv_gamma <- function(position, design, target, sigma = 100L){
  if(is.data.frame(design)) design <- as.matrix(design)
  dispersion <- position[1]
  beta <- position[-1]
  lin_predict <- design %*% beta
  comp_disp <- sum(target * exp(-lin_predict)) + sum(lin_predict)
  comp1 <- dispersion * t(design) %*% (target * exp(-lin_predict) - 1)
  c(comp_disp, comp1) - position/sigma
}


log_linear <- function(position, sigma = 10L){
  sum((-2*sigma^(-2))*(target - design %*% position)^2)
}

partial_deriv_lin <- function(position,sigma =10L){
  -sigma^(-2) * t(design) %*% (-target + design %*% position)
}
