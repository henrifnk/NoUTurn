build_efficient_tree_simply <- function(position_momentum, slice, direction, tree_depth,
                                 iter, is_log, stepsize, design = NULL, target = NULL, deltamax = 1e50,
                                 run = 1, count = 0, acceptance = 0, nacceptance = 0) {
  if(tree_depth == 0L) {# Basecase - take one leapfrogstep into direction
    build_leaf(position_momentum, slice, direction, tree_depth, stepsize, deltamax,
               design = design, target = target, iter = iter, is_log = is_log)
  } else {# Recursion- build left/right subtrees
    state <- build_efficient_tree(position_momentum, slice, direction, tree_depth - 1L,
                                  stepsize = stepsize, deltamax = deltamax, design = design,
                                  target = target, is_log = is_log, iter = iter,
                                  run = run, count = count,
                                  acceptance = acceptance, nacceptance = nacceptance)
    if(state$run){
      if(direction == -1L) {
        state1 <- build_efficient_tree(state$leftmost, slice, direction, tree_depth - 1L,
                                       stepsize = stepsize, deltamax = deltamax, design = design,
                                       target = target, is_log = is_log, iter = iter,
                                       run = run, count = count,
                                       acceptance = acceptance, nacceptance = nacceptance)
        state$leftmost <- state1$leftmost
        position_momentum <- state1$valid_state
      } else {
        state1 <- build_efficient_tree(state$rightmost, slice, direction, tree_depth - 1L,
                                       stepsize = stepsize, deltamax = deltamax, design = design,
                                       target = target, is_log = is_log, iter = iter,
                                       run = run, count = count,
                                       acceptance = acceptance, nacceptance = nacceptance)
        state$rightmost <- state1$rightmost
        position_momentum <- state1$valid_state
      }
      rate <- state1$count / (state1$count + state$count)
      if(!is.na(rate) && runif(1) <= rate) state$valid_state <- position_momentum
      # if any state is 0 Stopping criteria is triggered for trajectory iterations
      run <- run * state1$run * is_U_turn(state)
      count <- count + state$count + state1$count
      acceptance <- acceptance + state$acceptance$acceptance + state1$acceptance$acceptance
      nacceptance <- nacceptance + state$acceptance$nacceptance + state1$acceptance$nacceptance
      state$run <- run
      state$count <- count
      state$acceptance$acceptance <- acceptance
      state$acceptance$nacceptance <- nacceptance
    }

    return(state)
  }
}
# ______________________________________________________________________________
#' Build One Tree
#'
#' Builds one tree doing a leapfrogstep in one direction
#'
#' @inheritParams build_tree
build_leaf_simply <- function(position_momentum, slice, direction, tree_depth, stepsize, deltamax,
                       design = NULL, target = NULL, iter = NULL, is_log) {
  args <- if(is.null(design) && is.null(target)) {
    list(position_momentum$position)
  } else list(position_momentum$position, design, target)
  gradient_step <- do.call(gradient, args)
  step <- leapfrog(position_momentum$position, position_momentum$momentum,
                   stepsize = (direction * stepsize), gradient = gradient_step)
  proposal_state <- do.call(initialize_state, c(step, "efficient" = !is.null(iter)))
  proposal_density <- do.call(joint_probability, c(step, list(design), list(target), "is_log" = is_log))
  valid <- slice <= exp(proposal_density)
  if(!is.null(iter)) {
    proposal_state$count <- as.numeric(valid)
    quotient_lik <- exp(proposal_density - do.call(joint_probability, c(iter, list(design), list(target), "is_log" = is_log)))
    proposal_state$acceptance$acceptance <- min(1, quotient_lik)
    proposal_state$valid_state <- step
  }
  proposal_state$run <- 0 + (exp(proposal_density + deltamax) > slice)
  proposal_state
}
