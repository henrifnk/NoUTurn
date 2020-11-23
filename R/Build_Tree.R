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
build_tree <- function(position_momentum, slice, direction,
                       tree_depth, stepsize, deltamax = 1000L) {
  if(tree_depth == 0L) {# Basecase - take one leapfrogstep into direction
    build_leaf(position_momentum, slice, direction, tree_depth, stepsize, deltamax)
    } else {# Recursion- build left/right subtrees
      state <- build_tree(position_momentum, slice, direction, tree_depth - 1L, stepsize)
      if(direction == -1L) {
        state1 <- build_tree(state$leftmost, slice, direction, tree_depth - 1L, stepsize)
        position_momentum  <- state$leftmost <- state1$leftmost
      } else {
        state1 <- build_tree(state$rightmost, slice, direction, tree_depth - 1L, stepsize)
        position_momentum  <- state$rightmost <- state1$rightmost
      }
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
build_leaf <- function(position_momentum, slice, direction, tree_depth, stepsize, deltamax) {
  step <- do.call(leapfrog, c(position_momentum, "stepsize" = (direction * stepsize)))
  proposal_state <- do.call(initialize_state, step)
  proposal_density <- do.call(joint_probability, step)
  if(slice <= proposal_density) proposal_state$valid_state <- step
  proposal_state$run <- 0 + (proposal_density > (log(slice) - deltamax))
  proposal_state
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
  state2$valid_state$position <- rbind(state1$valid_state$position, state2$valid_state$position)
  state2$valid_state$momentum <- rbind(state1$valid_state$momentum, state2$valid_state$momentum)
  state2$valid_state
}
