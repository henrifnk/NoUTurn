#' Simplified naive U-Turn-Sampler Showcase
#'
#' @param position_init initial value for theta
#' @param stepsize size \epsilon for leapfrog steps
#' @param iteration number of iterations to sample from
#' @param design a design matrix, if posterior is a model
#' @param target a target feature vector if posterior is a model
#' @param is_log if the likelihood is already logged
#'
sample_naiveNUTS <- function(position_init, stepsize, iteration, design = NULL, target = NULL, is_log = TRUE) {
  ###### Preparation ############################################
  # initial position will be our starting value for the algorithm
  # positions will be a list containing all of our sampled positions
  position <- position_init
  positions <- data.frame(matrix(ncol = length(position_init), nrow = iteration))
  ###### NUTS Iteration ##########################################
  # For `iteration` steps we do:
  # we sample a momentum, combined with the position this will be
  # our proposal for Leapfrog Steps,
  # we initialize `tree_depth` as identifier for the depth and doubling
  # we form a slice & sample uniformly from 0 to the joint probability at current
  # position and momentum
  # we also set up a state that contains information about valid_states
  # sampled from Leapfrogsteps, our leftmost and rightmost state of
  # our trajectory and 'run' a criterion that tells us if we made the U-turn
  # during our last doubling
  for(m in 1L:(iteration)) {
    momentum <- rnorm(length(position))
    slice <- runif(1L, max = exp(joint_probability(position, momentum, design, target)))
    state <- initialize_state(position, momentum)
    tree_depth <- 0L
    ###### Build Tree ##########################################
    # While our Leapfrogsteps didn't made a U-Turn we do:
    # 1) sample a direction to integrate Leapfrogsteps
    #   a) foreward in time (+1)
    #   b) backward in time (-1)
    # 2) build a tree in the previously chosen direction
    # 3) update our right-/leftmost lepfrog node
    # 4) Check if U-Turn was made between
    #   a) sub trees performed in build tree
    #   b) the whole tree
    while(state$run) {
      # 1 means forward, -1 means backward doubling
      direction <- sample(c(-1L, 1L), 1L)
      state_proposal <- if(direction == -1) {
        build_tree(state$leftmost, slice, direction, tree_depth, stepsize = stepsize,
                   design = design, target = target, is_log = is_log)
      } else{
        build_tree(state$rightmost, slice, direction, tree_depth, stepsize = stepsize,
                   design = design, target = target, is_log = is_log)
      }
      if(state_proposal$run) {
        if(direction == -1) {
          state$leftmost <- state_proposal$leftmost
          state$valid_state <- unite_valid_states(state_proposal, state)
        } else{
          state$rightmost <- state_proposal$rightmost
          state$valid_state <- unite_valid_states(state, state_proposal)
        }
      }
      state$run <- state_proposal$run * is_U_turn(state = state)
      tree_depth <- tree_depth + 1
    }
    # Use the Transition kernel:
    # Sample uniformly one of the valid states drawn through Leapfrog steps
    if(is.matrix(state$valid_state$position) ||
       is.data.frame(state$valid_state$position)) {
      row_id <- sample(seq_len(nrow(state$valid_state$position)), 1L)
      positions[m, ] <- state$valid_state$position[row_id, ]
      position = as.numeric(positions[m, ])
    } else {
      positions[m, ] <- position
    }
  }
  positions
}
