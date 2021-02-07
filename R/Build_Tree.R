build_leaf <- function(position_momentum, direction, stepsize, slice) {
  step <- leapfrog(position_momentum$position, position_momentum$momentum,
                   stepsize = (direction * stepsize))
  state <- initialize_states(step$position, step$momentum)
  dens <- joint_log_density(step$position, step$momentum)
  if(slice <= exp(dens)) {
    state$valid_state <- step
    state$count = 1L
  }
  state$run <- (dens - log(slice)) >= -1e3
  if(is.na(state$run)) {
    warning("Deltamax induced NaN, moving in low posterior regions?")
    state$run = 1L
  }
  state
}

build_leaf_da <- function(position_momentum, direction, stepsize, slice, dens_0) {
  step <- leapfrog(position_momentum$position, position_momentum$momentum,
                   stepsize = (direction * stepsize))
  state <- initialize_states(step$position, step$momentum)
  dens <- joint_log_density(step$position, step$momentum)
  if(slice <= exp(dens)) {
    state$valid_state <- step
    state$count = 1L
  }
  state$run <- (dens - log(slice)) >= -1e3
  state$acceptance <- min(1, exp(dens - dens_0))
  if(is.na(state$run)) {
    warning("Deltamax induced NaN, moving in low posterior regions?")
    state$run = 1L
  }
  state
}


build_tree <- function(position_momentum, direction, tree_depth, stepsize, slice){
  if(tree_depth == 0L) {
    build_leaf(position_momentum, direction, stepsize, slice)
  } else {
    states <- build_tree(position_momentum, direction, tree_depth - 1L,
                         stepsize, slice)
    if(states$run) { #<<
      if(direction == -1L) {
        states_prop <- build_tree(states$leftmost, direction,
                                  tree_depth - 1L, stepsize, slice)
        position_momentum  <- states$leftmost <- states_prop$leftmost
      } else {
        states_prop <- build_tree(states$rightmost, direction,
                                  tree_depth - 1L, stepsize, slice)
        position_momentum  <- states$rightmost <- states_prop$rightmost
      }
      tree_ratio <- states_prop$count / (states_prop$count + states$count)#<<
      if(is.na(tree_ratio)) tree_ratio <- 0
      if(rbinom(1, 1, tree_ratio)) {
        states$valid_state <- states_prop$valid_state #<<
      }
      states$count <- states_prop$count + states$count #<<
      states$run <- states_prop$run || is_U_turn(states)
    }
    return(states)
  }
}

build_tree_da <- function(position_momentum, direction, tree_depth,
                          stepsize, slice, dens_0){
  if(tree_depth == 0L) {
    build_leaf_da(position_momentum, direction, stepsize, slice, dens_0)
  } else {
    states <- build_tree_da(position_momentum, direction, tree_depth - 1L,
                            stepsize, slice, dens_0)
    if(states$run) { #<<
      if(direction == -1L) {
        states_prop <- build_tree_da(position_momentum, direction,
                                     tree_depth - 1L, stepsize, slice, dens_0)
        position_momentum  <- states$leftmost <- states_prop$leftmost
      } else {
        states_prop <- build_tree_da(position_momentum, direction,
                                     tree_depth - 1L, stepsize, slice, dens_0)
        position_momentum  <- states$rightmost <- states_prop$rightmost
      }
      tree_ratio <- states_prop$count / (states_prop$count + states$count)
      if(is.na(tree_ratio)) tree_ratio <- 0
      if(rbinom(1, 1, tree_ratio)) {
        states$valid_state <- states_prop$valid_state #<<
      }
      states$steps <- states$steps + states_prop$steps
      states$acceptance <- states$acceptance + states_prop$acceptance
      states$count <- states_prop$count + states$count #<<
      states$run <- states_prop$run || is_U_turn(states)
    }
    return(states)
  }
}
