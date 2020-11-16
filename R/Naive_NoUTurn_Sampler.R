#' Naive U-Turn-Sampler
#'
#' @param theta_init initial value for theta
#' @param joint_prob log of joint density probability
#' @param iteration number of iterations to sample from
#' @param seed set a random seed to guarantee reproducability
#' @example posterior_density = mvtnorm::dmvnorm
#' @return
#'
#TODO: calculate posterior from data up to normalization constant
naive_nouturn_sampler <- function(position_init, stepsize, iteration, seed = 123L) {
  position <- list(position_init)
  set.seed(seed)
  for(m in 1:iteration) {
    position_leftm <- position_rightm <- position[[m]]
    momentum <- momentum_leftm <- momentum_rightm <- rnorm(length(position[[m]]))
    slice <- runif(1L, max = joint_probability(position[[m]], momentum))
    # 1 means True, 0 means False
    run <- 1L
    tree_depth = 0L
    valid_states <- list("position" = position[[m]], "momtentum" = momentum)
    while(run) {
      # 1 means forward, -1 means backward doubling
      direction <- sample(c(-1, 1), 1L)



   }
  }
}

joint_probability <- function(position, momentum) {
  exp(log(posterior_density(position)) - 0.5 * t(momentum) %*% momentum)
}
