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
build_trees <- function(state, deltamax = 1000) {
  if(state$tree_depth == 0) {
    state <- build_1tree(state, deltamax) # Basecase - take one leapfrogstep into direction
    }else{
      # Recursion- build left/right subtrees
      state$tree_depth = state$tree_depth - 1
      state1 <- build_trees(state, deltamax)
      state2 <- build_trees(state1, deltamax)
      # if any state is 0 Stopping criteria is triggered for trajectory iterations
      state2$run <- state1$run * state2$run * is_U_turn(state2, "state_leftmost") * is_U_turn(state2, "state_rightmost")
      state2$valid_state$position <- rbind(state2$valid_state$position, state1$valid_state$position)
      state2$valid_state$momentum <- rbind(state2$valid_state$momentum, state1$valid_state$momentum)
      state2
  }
}

#' Build One Tree
#'
#' Builds one tree doing a leapfrogstep in one direction
#'
#' @inheritParams build_trees
build_1tree <- function(state, deltamax) {
  proposal_state <- if(state$direction == -1) {
    state$state_leftmost <- do.call(leapfrog, c(state$state_leftmost,
                                                "stepsize" = (state$direction * state$stepsize)))
  } else {
    state$state_rightmost <-  do.call(leapfrog, c(state$state_ritghmost,
                                                  "stepsize" = (state$direction * state$stepsize)))
  }
  proposal_density <- do.call(joint_probability, proposal_state)
  if(slice <= proposal_density) state$valid_state <- proposal_state
  state$run <- 0 + (proposal_density > (log(slice) - deltamax))
  state
}

#' Is a U Turn made?
#'
#' Investigates if trajectory makes a U-Turn (if >= 0)
#' @inheritParams build_trees
is_U_turn <- function(state, direction) {
  momentum <- state[[direction]]$momentum
  position_distance <- state$state_ritghmost$position -state$state_leftmost$position

  as.numeric(t(momentum) %*% position_distance) >= 0
}

