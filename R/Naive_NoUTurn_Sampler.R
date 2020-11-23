#' Naive U-Turn-Sampler
#'
#' @param theta_init initial value for theta
#' @param joint_prob log of joint density probability
#' @param iteration number of iterations to sample from
#' @param stepsize size for leapfrog steps
#' @param seed set a random seed to guarantee reproducability
#' @example posterior_density = mvtnorm::dmvnorm
#' @return
#'
#TODO: calculate posterior from data up to normalization constant
naive_nouturn_sampler <- function(position_init, stepsize = NULL, iteration, seed = 123L, plot = FALSE) {
  set.seed(seed)
  if(is.null(stepsize)) stepsize <- find_initial_stepsize(position_init)
  position <- position_init
  positions <- vector("list", iteration)
  for(m in 1L:(iteration)) {
    momentum <- rnorm(length(position))
    slice <- runif(1L, max = joint_probability(position, momentum))
    # Initialize state to call on Build Tree
    state <- initialize_state(position, momentum)
    tree_depth <- 0L
    setup <- if(plot) setup_trajectory(position) else NULL
    while(state$run) {
      # 1 means forward, -1 means backward doubling
      direction <- sample(c(-1L, 1L), 1L)
      state_proposal <- if(direction == -1) {
        build_tree(state$leftmost, slice, direction, tree_depth, stepsize)
      } else{
        build_tree(state$rightmost, slice, direction, tree_depth, stepsize)
      }
      if(!is.null(setup)) setup <- plot_proposal(state_proposal$valid_state$position, setup)
      if(state_proposal$run) {
        if(direction == -1) {
          if(is.matrix(state_proposal$valid_state$position)) {
            # Sort backward doubled trajectory reversely
            state_proposal$valid_state <- lapply(state_proposal$valid_state, function(x) apply(x, 2, rev))
          }
          state$leftmost <- state_proposal$leftmost
          state$valid_state <- unite_valid_states(state_proposal, state)
        } else{
          state$valid_state <- unite_valid_states(state, state_proposal)
          state$rightmost <- state_proposal$rightmost
        }
      }
      state$run <- state_proposal$run * is_U_turn(state = state)
      tree_depth <- tree_depth + 1
      state$valid_state
      print(tree_depth)
    }
    row_id <- sample(seq_len(nrow(state$valid_state$position)), 1L)
    positions[[m]] = state$valid_state$position[row_id, ]
  }
}
#_______________________________________________________________________________
#'  Joint Probability of position and momentum
#'  @inheritParams leapfrog
#'
joint_probability <- function(position, momentum) {
  as.numeric(exp(log(posterior_density(position)) - 0.5 * t(momentum) %*% momentum))
}
#_______________________________________________________________________________
#' Initialize state
#'
#' Helperfunction to intialize state
#' @inheritParams leapfrog
#' @param slice slice sample drawn in U-Turn-Sampler
#' @inheritParams naive_nouturn_sampler
initialize_state <- function(position, momentum, run = 1L) {
  list(
    "valid_state" = structure(vector("list", length = 2L), names= c("position", "momentum")),
    "rightmost"= list("position" = position, "momentum" = momentum),
    "leftmost" = list("position" = position, "momentum" = momentum),  "run" = run
  )
}
