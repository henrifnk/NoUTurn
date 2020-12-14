#' Hamiltonian (hybrid) Monte Carlo
#'
#' @param position_init initial sample
#' @param stepsize length of Leapfrog stepsize
#' @param leapfrogsteps interger; amount of leapfrogsteps
#' @param iterations of HMC algorithm, should be long enough to draw samples from
#'
hamiltonianMC <- function(position_init, stepsize, leapfrogsteps,
                          iterations) {
  ###### Preparation ############################################
  # initial position will be our starting value for the algorithm
  # positions will be a list containing all of our sampled positions
  position <- position_init
  positions <- vector(mode = "list", length = iterations)
  ###### HMC Iteration ###########################################
  # For `iteration` steps we do:
  # we sample a momentum, combined with the position this will be
  # our proposal for Leapfrog Steps
  for(iter in seq_len(iterations)) {
    momentum <- rnorm(length(position_init))
    proposal <- list("position" = position, "momentum" = momentum)
    ###### Leapfrog Iteration #####################################
    # For `leapfrogsteps` steps we do:
    # We calculate the partial gradient of our log-posterior at the
    # drawn position & recalculate the next proposal by one more
    # leapfrogstep
    for(step in seq_len(leapfrogsteps)){
      gradient <- gradient(proposal$position)
      proposal <- do.call(leapfrog, c(proposal, stepsize, gradient))
    }
  # Like in Metropolis-Hastings, we compare the (unnormalized) posterior
  # density at our proposal with the density at our previously drawn
  # position in the ratio of those two estimates, which will form the
  # acceptance probability of our drawn proposal
  # If excepted we add the proposal to our output positions variable
    acceptance <- (joint_probability(proposal$position, proposal$momentum) /
                     joint_probability(position, momentum))
    acceptance <- min(1, acceptance)
    positions[[iter]] <- position <- if(runif(1) < acceptance) proposal$position else position
  }
  return(positions)
}
