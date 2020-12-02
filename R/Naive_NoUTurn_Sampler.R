#' Naive U-Turn-Sampler
#'
#' @param position_init initial value for theta
#' @param stepsize size \epsilon for leapfrog steps
#' @param iteration number of iterations to sample from
#' @param seed set a random seed to guarantee reproducability
#' @param plot if plot of trajectory should be done.
#' This only makes sense for presentation matters as plotting takes a long time.
#' @param design a design matrix, if posterior is a model
#' @param target a target feature vector if posterior is a model
#' @param is_log if the likelyhood is already logged
#' @example posterior_density = mvtnorm::dmvnorm
#' posterior_density = sigmoid_posterior
#' @return
#'
naive_nouturn_sampler <- function(position_init, stepsize, iteration, seed = 123L, plot = FALSE,
                                  design = NULL, target = NULL, is_log = TRUE) {
  set.seed(seed)
  position <- position_init
  positions <- data.frame(matrix(ncol = length(position_init), nrow = iteration))
  tree_depths <- vector("numeric", length = iteration)
  for(m in 1L:(iteration)) {
    momentum <- rnorm(length(position))
    slice <- runif(1L, max = exp(joint_probability(position, momentum, design, target, is_log = is_log)))
    # Initialize state to call on Build Tree
    state <- initialize_state(position, momentum)
    tree_depth <- 0L
    setup <- if(plot) setup_trajectory(position, positions) else NULL
    while(state$run) {
      # 1 means forward, -1 means backward doubling
      direction <- sample(c(-1L, 1L), 1L)
      state_proposal <- if(direction == -1) {
        build_tree(state$leftmost, slice, direction, tree_depth, stepsize, design, target)
      } else{
        build_tree(state$rightmost, slice, direction, tree_depth, stepsize, design, target)
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
          state$rightmost <- state_proposal$rightmost
          state$valid_state <- unite_valid_states(state, state_proposal)
        }
      }
      state$run <- state_proposal$run * is_U_turn(state = state)
      tree_depth <- tree_depth + 1
    }
    tree_depths[m] <- tree_depth
    if(is.matrix(state$valid_state$position) || is.data.frame(state$valid_state$position)) {
      row_id <- sample(seq_len(nrow(state$valid_state$position)), 1L)
      positions[m, ] <- state$valid_state$position[row_id, ]
      position = positions[m, ]
    } else {
      positions[m, ] <- position
    }
  }
  cbind(positions, "tree_depth" = tree_depths)
}
#_______________________________________________________________________________
#'  Joint Probability of position and momentum
#'  @inheritParams leapfrog
#'
joint_probability <- function(position, momentum, design, target, is_log = TRUE) {
  args <- if(is.null(design) && is.null(target)) list(position) else list(position, design, target)
  dens_estimate <- ifelse(is_log, do.call(posterior_density, args),
                          log(do.call(posterior_density, args))
                          )
  as.numeric(dens_estimate - 0.5 * t(momentum) %*% momentum)
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
