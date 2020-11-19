#' Build Trees doubles the postitions and momentums and such the trajectory length
#'
#' Repatedly called by NUTS until run == FALSE, that is when a U-Turn is performed.
#'
#' @param position_leftm leftmost tajectory position state
#' @param position_rightm rightmost tajectory position state
#' @param momentum_leftm leftmost tajectory momentum state
#' @param momentum_rightm rightmost tajectory momentum state
#' @param slice drawn slice sample to validate position-/momentum steps
#' @param direction determines whether trajectory should explore foreward or backwards
#' @param tree_depth balanced treedepth with momentum and position states in the tree nodes
#' @param stepsize parsed spesize \epsilon to Leapfrog step function
#'
#' @return
#'
build_tree <- function(state, tree_depth, deltamax = 1000L) {
  if(tree_depth == 0L) {
    build_leaf(state, deltamax) # Basecase - take one leapfrogstep into direction
    } else {
      # Recursion- build left/right subtrees
      state1 <- build_tree(state, tree_depth-1L, deltamax)
      state <- build_tree(state1, tree_depth-1L, deltamax)
      # if any state is 0 Stopping criteria is triggered for trajectory iterations
      state$run <- state1$run * state$run * is_U_turn(state)
      state$valid_state <- unite_valid_states(state1, state)
      state
  }
}
# ______________________________________________________________________________
#' Build One Tree
#'
#' Builds one tree doing a leapfrogstep in one direction
#'
#' @inheritParams build_tree
build_leaf <- function(state, deltamax) {
  proposal_state <- if(state$direction == -1) {
    state$leftmost <- do.call(leapfrog, c(state$leftmost,
                                          "stepsize" = (state$direction * state$stepsize)))
  } else {
    state$rightmost <-  do.call(leapfrog, c(state$rightmost,
                                            "stepsize" = (state$direction * state$stepsize)))
  }
  proposal_density <- do.call(joint_probability, proposal_state)
  if(state$slice <= proposal_density) state$valid_state <- proposal_state
  state$run <- 0 + (proposal_density > (log(state$slice) - deltamax))
  state
}
# ______________________________________________________________________________
#' Is a U Turn made?
#'
#' Investigates if trajectory makes a U-Turn (if >= 0)
#' @inheritParams build_tree
is_U_turn <- function(state, direction) {
  momentum_l <- state$leftmost$momentum
  momentum_r <- state$rightmost$momentum
  position_distance <- state$rightmost$position - state$leftmost$position
  left <- as.numeric((t(momentum_l) %*% position_distance)) >= 0
  right <- as.numeric((t(momentum_r) %*% position_distance)) >= 0
  left * right
}
# ______________________________________________________________________________
#' Unite states
#'
#' unites 2 valid state sets in second state set
unite_valid_states <- function(state1, state2) {
  state2$valid_state$position <- unique(rbind(state1$valid_state$position, state2$valid_state$position))
  state2$valid_state$momentum <- unique(rbind(state1$valid_state$momentum, state2$valid_state$momentum))
  state2$valid_state
}
