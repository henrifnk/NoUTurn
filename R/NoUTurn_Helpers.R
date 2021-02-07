#_______________________________________________________________________________
#' Initialize state
#'
#' Helperfunction to intialize state
#' @inheritParams leapfrog
#' @param run is u-turn made?
initialize_states <- function(position, momentum, run = 1L) {
  list(
    "valid_state" = structure(list(matrix(,nrow = 0, ncol = length(position)),
                                   matrix(,nrow = 0, ncol = length(position))),
                              names= c("position", "momentum")),
    "rightmost"= list("position" = position, "momentum" = momentum),
    "leftmost" = list("position" = position, "momentum" = momentum),
    "run" = run, "count" = 0, "acceptance" = 0, "steps" = 1
  )
}

# ______________________________________________________________________________
#' Is a U Turn made?
#'
#' Investigates if trajectory makes a U-Turn (if >= 0)
#' @param state state variable containing rightmost and leftmost position/doubling
is_U_turn <- function(state) {
  momentum_l <- state$leftmost$momentum
  momentum_r <- state$rightmost$momentum
  if(anyNA(c(state$rightmost$momentum, state$leftmost$momentum))){
    warning("momentum contains NA, trajectory aborted")
    return(FALSE)
  }
  position_distance <- state$rightmost$position - state$leftmost$position
  left <- sum(momentum_l * as.numeric(position_distance)) >= 0
  right <- sum(momentum_r * as.numeric(position_distance)) >= 0
  left * right
}
#_______________________________________________________________________________
#' Joint logarithmic density
#'
#' joint log density of position and momentum
#' @inheritParams leapfrog
joint_log_density <- function(position, momentum) {
  log_dens_estimate <- log_posterior_density(position)
  joint_dens <- log_dens_estimate - 0.5 * sum(momentum * momentum)
  if(is.na(joint_dens)) {
    warning("joint density induced NAs, density is set to -Inf")
    return(-Inf)
  }
  joint_dens
}
