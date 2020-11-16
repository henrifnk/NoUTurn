#' Build Tree doubles the postitions and momentums and such the trajectory length
#'
#' Repatedly called by NUTS until run == FALSE, that is when a U-Turn is performed.
#'
#' @param position_leftm leftmost tajectory position state
#' @param position_rightm rightmost tajectory position state
#' @param momentum_leftm leftmost tajectory momentum state
#' @param slice drawn slice sample to validate position-/momentum steps
#' @param direction determines whether trajectory should explore foreward or backwards
#' @param tree_depth balanced treedepth with momentum and position states in the tree nodes
#' @param stepsize parsed spesize \epsilon to Leapfrog step function
#'
#' @return
#'
build_tree <- function(position_leftm, position_rightm, momentum_leftm, momentum_rightm,
                       slice, direction, tree_depth, stepsize, valid_states) {
  mom <- if(direction == -1) momentum_leftm else momentum_rightm
  pos <- if(direction == -1) position_leftm else position_rightm

  for(i in 1: 2^tree_depth) {
    # Base Case: take one Leapfrogstep into sampled direction
    position_mometum <- leapfrog(position = pos, momentum = mom, stepsize = direction * stepsize)
    if(joint_probability(position_mometum$position, position_mometum$momentum) > slice) {
      valid_states <- list(valid_states, position_mometum)
    }
  }


}
