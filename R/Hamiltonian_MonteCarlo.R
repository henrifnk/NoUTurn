#' Hamiltonian (hybrid) Monte Carlo
#'
#' In computational physics and statistics, the Hamiltonian Monte Carlo algorithm (also known as hybrid Monte Carlo),
#' is a Markov chain Monte Carlo method for obtaining a sequence of random samples which converge to being
#' distributed according to a target probability distribution for which direct sampling is difficult.
#' This sequence can be used to estimate integrals with respect to the target distribution (expected values).
#'
#' @inheritParams naive_nouturn_sampler
#' @param leapfrogsteps interger; amount of leapfrogsteps
#'
hamiltonianMC <- function(position_init, stepsize, leapfrogsteps, iterations) {
  position <- position_init
  for(iter in seq_len(iterations)) {
    momentum <- rnorm(length(position_init))
    gradient
    proposal <- list("position" = position, "momentum" = momentum)
    for(step in seq_len(leapfrogsteps)){
      gradient <- gradient(proposal$position)
      proposal <- do.call(lepfrog, c(proposal, stepsize, gradient))
    }
    acceptance <- joint_probability(proposal$position, proposal$momentum) / joint_probability(position, momentum)
    acceptance <- min(1, acceptance)
    position <- sample(c(proposal$position, position), prob = c(acceptance, 1 - acceptance))
  }
}
