#' Efficient No-UTurn
#'
#' Efficient No-U-turn Sampler, implemented by Gelman et al. 2013
#' it is possible to determine stepsize by dual averaging
#'
#' @inheritParams naive_nouturn_sampler
#' @param adaption iterations where  stepsize is adapted. Must be smaller than iterations.
#' @param target_accpentance targeted average acceptance rate of proposals
#'
sample_noUturn <- function(position_init, iteration, seed = 1234L, stepsize = NULL,
                           adaption, target_accpentance = 0.65, design = NULL, target = NULL, is_log = TRUE) {
  set.seed(seed)
  if(is.null(stepsize)) stepsize <- find_initial_stepsize(position_init)
  duala_params <- init_parmeters(stepsize) # This is made fromheuristics, params can be adapted in init_params
  if(is.data.frame(design)) design <- as.matrix(design)
  position <- position_init
  positions <- data.frame(matrix(ncol = length(position_init), nrow = iteration))
  tree_depths <- vector("numeric", length = iteration)
  for(m in seq_len(iteration)) {
    momentum <- rnorm(length(position))
    slice <- runif(1L, max = exp(joint_probability(position, momentum, design, target, is_log = is_log)))
    # Initialize state to call on Build Tree
    state <- initialize_state(position, momentum, efficient = TRUE)
    iter <- list("position" = position, "momentum" = momentum)
    tree_depth <- 0L
    while(state$run) {
      # 1 means forward, -1 means backward doubling
      direction <- sample(c(-1L, 1L), 1L)
      state_proposal <- if(direction == -1) {
        build_efficient_tree(state$leftmost, slice, direction, tree_depth, stepsize, design, target, iter = iter)
      } else{
        build_efficient_tree(state$rightmost, slice, direction, tree_depth, stepsize, design, target, iter = iter)
      }
      if(state_proposal$run) {
        rate <- min(state_proposal$count / (state_proposal$count + state$count), 1)
        if(runif(1)<= rate) state$valid_state <- state_proposal$valid_state
      }
      state$count <- state_proposal$count + state$count
      state$run <- state_proposal$run * is_U_turn(state = state)
      tree_depth <- tree_depth + 1
    }
    if(m <= adaption) {
      duala_params <- update_stepsize(duala_params, m, target_accpentance,
                                      state$acceptance$acceptance / state$acceptance$nacceptance)
      stepsize <- duala_params$stepsize[m + 1]
    }
    if(m == adaption) {
      stepsize <- duala_params$stepsize_weight[m + 1]
    }
    tree_depths[m] <- tree_depth - 2L
    if(!is.null(state$valid_state$position)) position <- as.numeric(state$valid_state$position)
    positions[m, ] <- position
  }
  structure(cbind(positions, "tree_depth" = tree_depths),
            "dual_averaging" = duala_params)
}

#' Initialize parameters for dual averaging
#'
#' @param level where stepsize is shrunken towards (\mhy)
#' @param stepsize_weight weighted stepsize from previous iteration
#' @param mcmc_behavoir behavior of algorithm at previous iterations
#' @param shrinkage controls the amount of shrinkage
#' @param stability controls stability at early iterations
#' @param adaption controls the speed of converges and so the influence of more recent weights
#'
init_parmeters <- function(stepsize, level = NULL, stepsize_weight = 1,
                           mcmc_behavior = 0, shrinkage = 0.05,
                           stability = 10, adaption = 0.75) {
  if(is.null(level)) level <- log(10 * stepsize)

  list(
    "stepsize" = stepsize,
    "level" = level, "stepsize_weight" = stepsize_weight,
    "mcmc_behavior" = mcmc_behavior, "shrinkage" = shrinkage,
    "stability" = stability, "adaption" = adaption
  )
}

#' update stepsize
#'
#' updates stpesize parameters by performing dual averaging
#'
#' @inheritParams sample_noUturn
#' @param dap dual_avereging_params parameters necessary for dual averaging
#' @param iter actual iteration cycle
#' @param average_acceptance average achieved acceptance in iteration m
#'
update_stepsize <- function(dap, iter, target_accpentance, average_acceptance) {
  weight <- (iter + stability)^(-1)
  dap$mcmc_behavior[iter + 1] = (1 - weight) * dap$mcmc_behavior[iter] + weight * (target_accpentance - average_acceptance)
  dap$stepsize[iter + 1] <- exp(level - sqrt(iter) / shrinkage * dap$mcmc_behavior[iter + 1])
  dap$stepsize_weight[iter + 1] <- exp(iter^(-adaption) * log(dap$stepsize[iter + 1]) + (1 - iter^(-adaption)) * dap$stepsize[iter])
  dap
}
