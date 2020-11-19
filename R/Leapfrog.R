#' Lepfrog takes one step in trajectory
#'
#' Repatedly called by NUTS until run == FALSE, that is when a U-Turn is performed.
#'
#' @param position tajectory position state
#' @param momentum tajectory momentum state
#' @param stepsize stepsize \epsilon to Leapfrog step function
#'
#' @return
#'
leapfrog <- function(position, momentum, stepsize) {
  gradient <- gradient_approx(position)
  momentum <- momentum + (stepsize / 2) * gradient
  position <- position + stepsize * momentum
  momentum <- momentum + (stepsize / 2) * gradient
  return(list("position" = position, "momentum" = momentum))
}
#_______________________________________________________________________________
#' Gradient Approx
#'
#' simple approximation scheme for gradient
#'
#' @param x x-value where derivation should be evaluated
#' @param delta absolute range to each side of x
#' @param n size of samples where function should be evaluated
#'
gradient_approx = function(x, delta = 1e-5, n = 3){
 x <- sapply(x, function(x) seq(from = x - delta, to = x + delta, length.out = max(2, n)))
 y <- apply(x, MARGIN = 1, function(x) posterior_density(x))
 mean(diff(y)/diff(x))
}
