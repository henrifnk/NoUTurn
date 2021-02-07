#' Hamiltonian (hybrid) Monte Carlo
#'
#' @param position_init initial sample
#' @param stepsize length of Leapfrog stepsize
#' @param leapfrogsteps interger; amount of leapfrogsteps
#' @param iterations of HMC algorithm, should be long enough to draw samples from
#'
hamiltonianMC <- function(position_init, stepsize, leapfrogsteps, iterations) {
  pb <- txtProgressBar(min = 0, max = iterations, style = 3)
  ###### Preparation ############################################
  # initial position will be our starting value for the algorithm
  # positions will be a list containing all of our sampled positions
  position <- position_init
  positions <- data.frame(matrix(ncol = length(position_init), nrow = iterations))
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
      proposal <- leapfrog(proposal$position, proposal$momentum, stepsize)
    }
  # Like in Metropolis-Hastings, we compare the (unnormalized) posterior
  # density at our proposal with the density at our previously drawn
  # position in the ratio of those two estimates, which will form the
  # acceptance probability of our drawn proposal
  # If excepted we add the proposal to our output positions variable
    acceptance <- exp(joint_log_density(proposal$position, proposal$momentum) -
                     joint_log_density(position, momentum))
    if(is.na(acceptance)) acceptance <- 1
    acceptance <- min(1, acceptance)
    position <- if(runif(1) < acceptance) proposal$position else position
    positions[iter, ] <- position
    setTxtProgressBar(pb, iter)
  }
  return(positions)
}
