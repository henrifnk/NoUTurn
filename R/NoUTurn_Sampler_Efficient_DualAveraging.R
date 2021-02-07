#' Sample No-u-Turn with dual averaging
#'
#' NUTS with dual averaging introduced by Hoffman and Gelman
#'
#' @inheritParams sample_NUT
#' @param warmup MCMC iterations that shall be used for warm up phase
#' @param target_accpentance target acceptance for dual averaging to estimate stepsize
#' @param max_depth maximal tree_depth for an MCMC iteration
#' @export
sample_NUT_da <- function(init_position, iteration, warmup, target_accpentance = 0.65,
                          max_depth = 13, seed = 123L){
  set.seed(seed)
  pb <- txtProgressBar(min = 0, max = iteration, style = 3)
  position <- init_position
  stepsize <- find_initial_stepsize(init_position) #<
  duala_param <- init_parmeters(stepsize)
  positions <- data.frame(matrix(ncol = length(init_position), nrow = iteration))
  for(iter in seq_len(iteration)){
    momentum <- rnorm(length(position))
    dens <- joint_log_density(position, momentum)
    if(is.na(dens)) {
      warning(paste("NUTS sampled NA in iteration", iter))
      dens = 1
    }
    slice <- runif(n = 1, min = 0, max = exp(dens))
    states <- initialize_states(position, momentum)
    tree_depth <- 0L
    while(states$run){
      direction <- sample(c(-1, 1), 1)
      if(direction == -1L) {
        states_prop <- build_tree_da(states$leftmost, direction, tree_depth,
                                     stepsize, slice, dens)
        states$leftmost <- states_prop$leftmost
      } else {
        states_prop <- build_tree_da(states$rightmost, direction, tree_depth,
                                     stepsize, slice, dens)
        states$rightmost <- states_prop$rightmost
      }
      if(states_prop$run) {
        states$run <- is_U_turn(states) * states_prop$run
        tree_ratio <- min(1, states_prop$count / states$count)
        if(is.na(tree_ratio)) tree_ratio <- 0
        if(rbinom(1, 1, tree_ratio)) {
          states$valid_state <- states_prop$valid_state
        }
      } else break
      states$count <- states_prop$acceptance + states_prop$count #<<
      if(tree_depth == max_depth) break
      tree_depth = tree_depth + 1
    }
    if(iter <= warmup)  {
      average_acceptance <- states_prop$acceptance / states_prop$step
      duala_param <- update_stepsize(duala_param, iter,
                                     target_accpentance, average_acceptance)
      stepsize <- duala_param$stepsize[iter + 1]
    }
    if(iter == (warmup + 1)) stepsize <- duala_param$stepsize_weight[iter]
    setTxtProgressBar(pb, iter)
    if(is.matrix(states$valid_state$position)) {
      positions[iter,] <- position
      next
    }
    position <- as.numeric(states$valid_state$position)
    positions[iter,] <- position
  }
  return(positions)
}
