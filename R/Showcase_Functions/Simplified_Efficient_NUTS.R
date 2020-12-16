#' Efficient No-UTurn
#'
#' Efficient No-U-turn Sampler, implemented by Gelman et al. 2013
#' it is possible to determine stepsize by dual averaging
#'
#' @inheritParams naive_nouturn_sampler
#' @param adaption iterations where  stepsize is adapted. Must be smaller than iterations.
#' @param target_accpentance targeted average acceptance rate of proposals
#'
sample_noUturn <- function(position_init, iteration,
                           adaption = 0.5 * iteration, target_accpentance = 0.65,
                           design = NULL, target = NULL, is_log = TRUE, max_tree_depth = 10L) {
  stepsize <- find_initial_stepsize(position_init, design, target, is_log)
  duala_params <- init_parmeters(stepsize)
  position <- position_init
  positions <- data.frame(matrix(ncol = length(position_init), nrow = iteration))
  for(m in seq_len(iteration)) {
    momentum <- rnorm(length(position))
    iter <- list("position" = position, "momentum" = momentum)
    slice <- runif(1L, max = exp(joint_probability(position, momentum, design,
                                                   target, is_log = is_log)))
    # Initialize state to call on Build Tree
    state <- initialize_state(position, momentum, efficient = TRUE)
    state$count <- 1
    tree_depth <- 0L
    while(state$run) {
      # 1 means forward, -1 means backward doubling
      direction <- sample(c(-1L, 1L), 1L)
      if(direction == -1) {
        state_proposal <- build_efficient_tree(state$leftmost, slice, direction,
                                               tree_depth, stepsize = stepsize, design = design,
                                               target = target, is_log = is_log, iter = iter)
        state$leftmost <- state_proposal$leftmost
      } else{
        state_proposal <- build_efficient_tree(state$rightmost, slice, direction,
                                               tree_depth, stepsize = stepsize, design = design,
                                               target = target, is_log = is_log, iter = iter)
        state$rightmost <- state_proposal$rightmost
      }
      if(state_proposal$run) {
        rate <- min(state_proposal$count / state$count, 1)
        if(rbinom(n = 1, size = 1, prob = rate)) state$valid_state <- state_proposal$valid_state
      }
      state$count <- state_proposal$count + state$count
      state$run <- state_proposal$run * is_U_turn(state = state)
      if(tree_depth > max_tree_depth){
        warning("NUTS: Reached max tree depth")
        break
      }
      tree_depth <- tree_depth + 1
    }
    if(m <= adaption) {
      duala_params <- update_stepsize(duala_params, m, target_accpentance,
                                      state_proposal$acceptance$acceptance / state_proposal$acceptance$nacceptance)
      stepsize <- duala_params$stepsize[m + 1]
    }
    if(m == adaption) {
      stepsize <- duala_params$stepsize_weight[m + 1]
    }
    tree_depths[m] <- tree_depth - 2L
    if(!is.null(state$valid_state$position)) position <- as.numeric(state$valid_state$position)
    positions[m, ] <- position
  }
  structure(positions, "tree_depth" = tree_depths,
            "dual_averaging" = duala_params)
}
