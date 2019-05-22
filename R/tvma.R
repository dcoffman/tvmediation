
tvma <- function(treatment, t.seq, mediator, outcome, t.est = t.seq, plot = TRUE, CI="boot", replicates = 1000) {
  # Calculate time-varying mediation effect and sd
  #
  # Args:
  #   treatment   -->   a vector with treatment values
  #   t.seq       -->   a vector of time points for each observation
  #   mediator    -->   matrix of mediator values in wide format
  #   outcome     -->   matrix of outcome outcomes in wide format
  #   t.est       -->   time points to make the estimation                              Defualt = t.seq
  #   plot        -->   TRUE or FALSE for plotting mediation effect                     Default = "FALSE".
  #   CI          -->   "asymptotic" or "boot" method of deriving confidence intervals. Default = "boot".
  #   replicates  -->   Number of replicates for bootstrapping confidence intervals.    Default = 1000.
  #
  # Returns:
  #   hat.alpha.1       -->   estiamted treatment effect on mediator
  #   hat.beta.2        -->   estiamted mediation effect on outcome
  #   mediation.effect  -->   time varying mediation effect
  #   sd                -->   estimated standard deviation of the mediation effect (only for asymptotic estimation)
  #   CI.upper          -->   Upper confidence intervals
  #   CI.lower          -->   Lower confidence intervals
  #
  ##

  deltat <- max(diff(t.seq)) / 2  # half the time between two measures
  N <- length(treatment)
  t.coeff <- NULL
  for (j in 2:length(t.seq)) {
    # create empty vector, store raw mediator.
    temp.mediator.j  <- NULL
    temp.mediator.j  <- cbind(mediator[j - 1, ], mediator[j, ])
    # Derive centered Mediators and Outcomes
    newMO.j.est <- newMediatorOutcome(treatment, temp.mediator.j, outcome[j - 1, ])
    # Estimate coefficients, then store them.
    coeff.est <- estCoeff(newMO.j.est)
    t.coeff <- cbind(t.coeff, coeff.est)  # store coeff estimates at t.seq
  }

  # EQUATIONS 4 & 5
  est.smooth <- smoothest(t.seq, t.coeff, t.est, deltat)

  # CALCULATE CONFIDENCE INTERVALS
  if (CI == "asymptotic") {
    results <- estAsymptoticCIs(treatment, t.seq, mediator, outcome, t.est, t.coeff, est.smooth, deltat)
    plot.label <- "Asymptotic CI"
  } else if (CI == "boot") {
    results <- estBootCIs(treatment, t.seq, mediator, outcome, t.est, deltat, replicates)
    plot.label <- "Bootstrap CI"
  }
  results <- c(list(hat.alpha.1 = est.smooth$hat.alpha.1, hat.beta.2 = est.smooth$hat.beta.2,
              mediation.effect = est.smooth$true.M), results)


  ## plot the time-varying mediation effect and the 95% confidence band
  if (plot == TRUE) {
    y.min <- min(results$CI.lower, na.rm = TRUE)
    y.max <- max(results$CI.upper, na.rm = TRUE)
    x.min <- floor(min(t.est))
    x.max <- ceiling(max(t.est))
    plot(t.est, results$mediation.effect, type = "l", col = "blue",
         main = "Time Varying Mediation",
         ylab = "Mediation Effect", ylim = c(y.min, y.max),
         xlab = "Time Points",      xlim = c(x.min, x.max), xaxt = 'n')
    axis(1, at = t.seq)
    lines(t.est, results$CI.upper, col = 2, lty = 2)
    lines(t.est, results$CI.lower, col = 2, lty = 2)
    legend("topright", inset = 0.02, legend = c("tvm", plot.label), col = c("blue", "red"), lty = 1:2, box.lty = 0)
  }

  ## Print results to screen and return them
  print("Time Varying Mediation Results:")
  print(results)
  return(results)

}
