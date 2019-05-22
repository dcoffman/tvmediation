
kern <- function(x, h, x.obs) {
  # Kernal function...
  #
  # Args:
  #   x     -->   t.est at time j or t-deltat
  #   h     -->   smoothed alpha or beta
  #   x.obs -->   a vector of time points at each obs[-1]
  #
  # Returns:
  #   kern  -->
  #
  ##

  del <- x - x.obs
  ndel <- del^2
  w <- ((1 / sqrt(2 * pi))) * (1 / h) * exp((-1 / 2) * ndel / h^2)
  return(w)

}
