#' Set Up A plot for trajectory
#'
setup_trajectory <- function(position, valid_positions, fun_dim1 = dnorm, fun_dim2 = dnorm, limit_dim1 = c(-3, 3), limit_dim2 = c(-3, 3)){
  library(ggplot2)
  library(gridExtra)
  theme_set(theme_classic() +  theme(axis.title.x = element_blank(),
                                     axis.title.y = element_blank(),
                                     axis.text.x = element_blank(),
                                     axis.text.y = element_blank(),
                                     axis.ticks = element_blank(),
                                     axis.line = element_blank())
            )
  density1 <- ggplot(data = data.frame(x = limit_dim1), aes(x)) +
    geom_area(stat = "function", fun = fun_dim1, fill = "blue", alpha =.5, args = list(mean = 0, sd = 1))
  density2 <- ggplot(data = data.frame(x = limit_dim2), aes(x)) +
    geom_area(stat = "function", fun = fun_dim2, fill = "blue", alpha =.5, args = list(mean = 0, sd = 1)) +
    coord_flip()

  blank_plot <- ggplot() + geom_blank(aes(1, 1)) +
    theme(
      plot.background = element_blank(),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank()
    )

  x <- seq(from = -3, to = 3, length.out = 20)
  grid <- expand.grid(x1 = x, y1 =x)
  grid <- cbind(grid, z1 = apply(grid, 1, posterior_density))

  valid_positions = na.omit(valid_positions)
  if(!is.data.frame(position)) position <- as.data.frame(t(position))
  position = setNames(position, c("a", "b"))

  scatter_plot <- ggplot(data = position, aes(a, b)) +
    stat_contour(alpha = 0.5, data = grid, aes(x = x1, y = y1, z = z1)) +
    geom_point(size = 1L) +
    xlim(limit_dim1) +
    ylim(limit_dim2)

  if(nrow(valid_positions) > 0L) {
    scatter_plot <- scatter_plot + geom_point(data = valid_positions, aes(X1, X2), colour = 'green')
  }
  list(density1, blank_plot, scatter_plot, density2)

}
#' plot trajectory induced by leapfrog
#'
#'
plot_proposal <- function(proposal, setup) {
  if(!is.null(proposal)) {
    if(is.null(dim(proposal))) {
      names(proposal) <- c("a", "b")
    } else colnames(proposal) <- c("a", "b")
    setup[[3]]$data <- rbind(setup[[3]]$data, proposal)
  }
  grid.arrange(
    setup[[1]], setup[[2]], setup[[3]], setup[[4]],
    ncol=2, nrow=2, widths=c(4, 1), heights=c(1, 4)
  )

  Sys.sleep(2)
  setup
}
