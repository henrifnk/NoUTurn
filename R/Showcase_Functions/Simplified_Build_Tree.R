#' Build Trees doubles the positions and momentum and such the trajectory length
#'
#' Repeatedly called by NUTS until run == FALSE, that is when a U-Turn is performed.
#'
#' @param position_momentum list, containing actual position and momentum
#' @param slice drawn slice sample to validate position-/momentum steps
#' @param direction determines whether trajectory should explore forwards or backwards
#' @param tree_depth balanced tree depth with momentum and position states in the tree nodes
#' @param is_log logical, indicator if joint probability is already logged
#' @param stepsize parsed stepsize \epsilon to Leapfrog step function
#' @param design design matrix of features for linear modell
#' @param target target variable to validate predictions
#' @param deltamax stop if the error in simulation increases
#' @param run independent argument that tells build tree if stopping criteria is met
#' @param valid_states all valid states visited so far
#'
build_tree_simply <- function(position_momentum, slice, direction, tree_depth,
                              is_log, stepsize, design = NULL, target = NULL,
                              deltamax = 1000L, run = 1,
                              valid_states = structure(vector("list", length = 2L))) {
  ###### Base Case ###########################################
  # if tree_depth is zero: Call build leaf,
  # take one leapfrog step into direction
  if(tree_depth == 0L) {
    build_leaf(position_momentum, slice, direction, tree_depth,
               stepsize = stepsize, deltamax = deltamax, design = design,
               target = target, is_log = is_log)
  } else {
    ###### Recursion ###########################################
    # Recursion- build left/right subtrees
    # Call build_tree at tree_depth - 1 twice
    # leapfrogsteps will get called 2^treedepth times
    state <- build_tree(position_momentum, slice, direction, tree_depth - 1L,
                        stepsize = stepsize, deltamax = deltamax, design = design,
                        target = target, is_log = is_log, run = run,
                        valid_states = valid_states)
    # After calling the second build_tree, assign depending on
    # direction left or rightmost node to state, which will be returned
    if(direction == -1L) {
      state1 <- build_tree(state$leftmost, slice, direction, tree_depth - 1L,
                           stepsize = stepsize, deltamax = deltamax, design = design,
                           target = target, is_log = is_log, run = run,
                           valid_states = valid_states)
      state$leftmost <- state1$leftmost
    } else {
      state1 <- build_tree(state$rightmost, slice, direction, tree_depth - 1L,
                           stepsize = stepsize, deltamax = deltamax, design = design,
                           target = target, is_log = is_log, run = run,
                           valid_states = valid_states)
      state$rightmost <- state1$rightmost
    }
    ###### Return State #######################################
    # Check if any stopping criteria was met during recursion
    # update valid states by unioning states from both build
    # tree calls
    run <- run * state1$run * state$run * is_U_turn(state)
    state$run <- run
    valid_states$position <- rbind(valid_states$position,
                                   state$valid_state$position,
                                   state1$valid_state$position)
    valid_states$momentum <- rbind(valid_states$momentum,
                                   state$valid_state$momentum,
                                   state1$valid_state$momentum)
    state$valid_state$position <- valid_states$position
    state$valid_state$momentum <- valid_states$momentum
    return(state)
  }
}

build_leaf_simply <- function(position_momentum, slice, direction, tree_depth,
                              stepsize, deltamax, design = NULL, target = NULL,
                              is_log) {
  ###### Call Gradient #####################################
  # We call the posteriors gradient at the current state of
  # our position momentum
  # In case of a GRM we need to specify design and target
  args <- if(is.null(design) && is.null(target)) {
    list(position_momentum$position)
  } else list(position_momentum$position, design, target)
  gradient_step <- do.call(gradient, args)
  ###### Leapfrog Step ####################################
  # From our actual position:
  # We do one Leapfrog step to the sampled direction
  # This step is largely determined by the gradient
  step <- leapfrog(position_momentum$position,
                   position_momentum$momentum,
                   stepsize = (direction * stepsize),
                   gradient = gradient_step)
  ###### Return State #######################################
  # We create a new state containing our leapfrog step as
  # rightmost, leftmost and current state
  # We calculate the likelihood of our (unnormalized)
  # joint probability at our current state
  # if this is more likely than the slice drawn from state of
  # the previous iteration, current state is assigned as valid state
  # Additionally we check if simulation error didn't increase too much
  proposal_state <- do.call(initialize_state, step)
  proposal_density <- do.call(joint_probability, c(step, design, target,
                                                   "is_log" = is_log))
  valid <- slice <= exp(proposal_density)
  if(valid) proposal_state$valid_state <- step
  proposal_state$run <- 0 + (proposal_density > log(slice) - deltamax)
  proposal_state
}
