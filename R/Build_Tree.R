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
#' @param stepsize parsed stepsize \epsilon to Leapfrog step function
#'
#' @return
#'
build_tree <- function(position_momentum, slice, direction, tree_depth, is_log,
                       stepsize, design = NULL, target = NULL, deltamax = 1000L,
                       run = 1, valid_states = structure(vector("list", length = 2L))) {
  if(tree_depth == 0L) {# Basecase - take one leapfrogstep into direction
    build_leaf(position_momentum, slice, direction, tree_depth,
               stepsize = stepsize, deltamax = deltamax, design = design,
               target = target, is_log = is_log)
    } else {# Recursion- build left/right subtrees
      state <- build_tree(position_momentum, slice, direction, tree_depth - 1L,
                          stepsize = stepsize, deltamax = deltamax, design = design,
                          target = target, is_log = is_log, run = run, valid_states = valid_states)
      if(direction == -1L) {
        state1 <- build_tree(state$leftmost, slice, direction, tree_depth - 1L,
                             stepsize = stepsize, deltamax = deltamax, design = design,
                             target = target, is_log = is_log, run = run, valid_states = valid_states)
        position_momentum  <- state$leftmost <- state1$leftmost
      } else {
        state1 <- build_tree(state$rightmost, slice, direction, tree_depth - 1L,
                             stepsize = stepsize, deltamax = deltamax, design = design,
                             target = target, is_log = is_log, run = run, valid_states = valid_states)
        position_momentum  <- state$rightmost <- state1$rightmost
      }
      # if any state is 0 Stopping criteria is triggered for trajectory iterations
      run <- run * state1$run * state$run * is_U_turn(state)
      state$run <- run
      valid_states$position <- rbind(valid_states$position, state$valid_state$position, state1$valid_state$position)
      valid_states$momentum <- rbind(valid_states$momentum, state$valid_state$momentum, state1$valid_state$momentum)
      state$valid_state$position <- valid_states$position
      state$valid_state$momentum <- valid_states$momentum
      return(state)
  }
}

build_efficient_tree <- function(position_momentum, slice, direction, tree_depth,
                                 iter, is_log, stepsize, design = NULL, target = NULL, deltamax = 1e50L,
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
build_leaf <- function(position_momentum, slice, direction, tree_depth, stepsize, deltamax,
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
  proposal_state$run <- if(slice > 0) 0 + (exp(proposal_density + deltamax) > slice) else 1
  proposal_state
}
# ______________________________________________________________________________
#' Is a U Turn made?
#'
#' Investigates if trajectory makes a U-Turn (if >= 0)
#' @inheritParams build_tree
is_U_turn <- function(state) {
  momentum_l <- state$leftmost$momentum
  momentum_r <- state$rightmost$momentum
  position_distance <- state$rightmost$position - state$leftmost$position
  left <- as.numeric((t(momentum_l) %*% as.numeric(position_distance))) >= 0
  right <- as.numeric((t(momentum_r) %*% as.numeric(position_distance))) >= 0
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
