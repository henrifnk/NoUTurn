#' Naive U-Turn-Sampler
#'
#' @param theta_init initial value for theta
#' @param joint_prob log of joint density probability
#' @param iteration number of iterations to sample from
#' @param seed set a random seed to guarantee reproducability
#'
#' @return
#'
naive_nouturn_sampler <- function(position_init, stepsize, joint_prob, iteration, seed = 123) {
  position <- list("init" = position_init)
  position_leftm <- position_rightm <- position[["init"]]
  set.seed(seed)
  for(m in 1:iteration) {
   momentum <- momentum_leftm <- momentum_rightm <- rnorm(length(position[["init"]]))
   run <- TRUE
   while(run) {
     #TODO: calculate posterior from data up to normalization constant
     log_posterior <- log(mvtnorm::dmvnorm(theta))
     # 1 means foreward, -1 means backward doubling
     direction <- sample(c(-1, 1), 1L)


   }
  }
}
