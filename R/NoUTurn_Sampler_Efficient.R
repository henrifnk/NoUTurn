#' Sample No-u-Turn
#'
#' efficient NUTS introduced by Hoffman and Gelman
#'
#' @inheritParams  sample_nNUT
#' @export
sample_NUT <- function(init_position, stepsize, iteration, seed = 123L){
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
    tree_depth <- 0L
    while(states$run){
      direction <- sample(c(-1, 1), 1)
      if(direction == -1L) {
        states_prop <- build_tree(states$leftmost, direction,
                                  tree_depth, stepsize, slice)
        states$leftmost <- states_prop$leftmost
      } else {
        states_prop <- build_tree(states$rightmost, direction,
                                  tree_depth, stepsize, slice)
        states$rightmost <- states_prop$rightmost
      }
      if(states_prop$run) {
        if(is.na(is_U_turn(states))) browser()
        states$run <- is_U_turn(states) * states_prop$run
        tree_ratio <-min(1, states_prop$count / states$count) #<<
        if(is.na(tree_ratio)) tree_ratio <- 0
        if(rbinom(1, 1, tree_ratio)) {
          states$valid_state <- states_prop$valid_state #<<
        }
      } else break
      states$count <- states_prop$count + states_prop$count #<<
      tree_depth = tree_depth + 1
    }
    setTxtProgressBar(pb, iter)
    if(is.matrix(states$valid_state$position)) {
      positions[iter,] <- position
      next
    }
    position <- as.numeric(states$valid_state$position)
    if(anyNA(position)) stop("position contains NAs")
    positions[iter,] <- position
  }
  return(positions)
}
