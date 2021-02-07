#_______________________________________________________________________________
#' Naive No U Turn Sampler
#'
#' Sample positions with naive u turn sampler
#'
#' @param init_position vector; intial position to start MCMC
#' @param stepsize numeric; size of Leapfrog steps
#' @param iteration integer; amount of MCMC iterations
#' @param seed integer; for reproducability
#' @export
sample_nNUT <- function(init_position, stepsize, iteration, seed = 123L){
  set.seed(seed)
  pb <- txtProgressBar(min = 0, max = iteration, style = 3)
  position <- init_position
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
    slice <-
      tree_depth <- 0L
    while(states$run){
      direction <- sample(c(-1, 1), 1)
      if(direction == -1L) {
        states_prop <- build_tree(states$leftmost, direction, tree_depth, stepsize)
        states$leftmost <- states_prop$leftmost
      } else {
        states_prop <- build_tree(states$rightmost, direction, tree_depth, stepsize)
        states$rightmost <- states_prop$rightmost
      }
      states$run <- is_U_turn(states) * states$run * states_prop$run
      if(states$run) {
        states$valid_state$position <- rbind(states$valid_state$position,
                                             states_prop$valid_state$position)
        states$valid_state$momentum <- rbind(states$valid_state$momentum,
                                             states_prop$valid_state$momentum)
      }
      tree_depth = tree_depth + 1
    }
    setTxtProgressBar(pb, iter)
    valid_steps <- nrow(states$valid_state$position)
    if(valid_steps == 0) {
      positions[iter,] <- position
      next
    }
    position <- states$valid_state$position[sample(seq_len(valid_steps), 1),]
    positions[iter,] <- position
  }
  return(positions)
}
