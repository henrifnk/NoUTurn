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
leapfrog <- function(position, momentum, stepsize, gradient) {
  momentum <- momentum + (stepsize / 2) * gradient
  position <- position + stepsize * momentum
  momentum <- momentum + (stepsize / 2) * gradient
  return(list("position" = position, "momentum" = momentum))
}
