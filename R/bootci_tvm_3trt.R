#' Bootstrap samples to estimate confidence intervals for continuous outcome and three 
#' treatment groups.
#' 
#' Internal function for estimating bootstrapped confidence intervals for the mediation 
#' effect of continuous outcome and three treatment groups when user argument 
#' \code{CI="boot"}.
#' 
#' @param boot.sample    number of replicates for bootstrapping confidence intervals. 
#'                       Default = 1000.
#' @param orig.data      a list of original data \code{T1}, \code{T2}, \code{mediator}, 
#'                       \code{outcome} and \code{t.seq}.
#' @param t.est          time points at which to make the estimation. Default = t.seq.
#' 
#' @return \item{plw1}{CI lower limit for estimated mediation effect of T1}
#' @return \item{pup1}{CI upper limit for estimated mediation effect of T1}
#' @return \item{plw2}{CI lower limit for estimated mediation effect of T2}
#' @return \item{pup2}{CI upper limit for estimated mediation effect of T2}
#' @return \item{orig.se1.all}{estimated standard errors for the mediation effect of T1}
#' @return \item{orig.se2.all}{estimated standard errors for the mediation effect of T2}
#' @return \item{orig.mediation1}{time varying mediation effect for T1}
#' @return \item{orig.mediation2}{time varying mediation effect for T2}
#' 
#' 

bootci_tvm_3trt <- function(boot.sample, orig.data, t.est)
{
  N <- length(orig.data$T1)
  x <- orig.data$x
  y <- orig.data$y
  T1 <- orig.data$T1
  T2 <- orig.data$T2
  t.seq <- orig.data$t.seq
  deltat <- max(diff(t.seq))/2
  result.org <- tvmcurve_3trt(T1, T2, t.seq, x, y, t.est)

  original.mediation1 <- (result.org$hat.mediation1)
  original.mediation2 <- (result.org$hat.mediation2)
  
  storage.boot1 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot2 <- matrix(0, nrow=boot.sample, ncol=length(t.est))
  
  start.time <- Sys.time()
  print("Beginning bootstrap for mediation effect CIs.")
  
  for(k1 in 1:boot.sample)
  {
    if (k1 < boot.sample) {
      cat(sprintf("Boostrapping iteration %03d", k1), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", k1))
    }
    index.sample <- sample(1:N, N, replace=TRUE)
    x.boot <- x[,index.sample]
    y.boot <- y[,index.sample]
    T1.boot <- T1[index.sample]
    T2.boot <- T2[index.sample]
    deltat <- max(diff(t.seq))/2
    
    result.boot <- tvmcurve_3trt(T1.boot, T2.boot, t.seq, x.boot, y.boot, t.est) 
    storage.boot1[k1,] <- result.boot$hat.mediation1
    storage.boot2[k1,] <- result.boot$hat.mediation2
  }
  
  orig.se1.all <- apply(storage.boot1, 2, sd)
  orig.se2.all <- apply(storage.boot2, 2, sd)
  
  t.length <- dim(storage.boot1)[2]
  lower1 <- rep(0,t.length)
  upper1 <- rep(0,t.length)
  lower2 <- rep(0,t.length)
  upper2 <- rep(0,t.length)
  
  for(k2 in 1:(dim(storage.boot1)[2])) # run through different time points 
  {
    temp1 <- storage.boot1[,k2]
    temp2 <- storage.boot2[,k2]
    
    orig.est1 <- original.mediation1[k2]
    orig.est2 <- original.mediation2[k2]
    
    ###percentile bootstrap method 
    upper1[k2] <- quantile(temp1,0.975)
    lower1[k2] <- quantile(temp1,0.025)
    
    upper2[k2] <- quantile(temp2,0.975)
    lower2[k2] <- quantile(temp2,0.025)
  }

  results <- list(plw1=lower1, pup1=upper1, plw2=lower2, pup2=upper2, 
                  orig.se1.all=orig.se1.all, orig.se2.all=orig.se2.all, 
                  orig.mediation1=original.mediation1, orig.mediation2=original.mediation2) 
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", 
                as.numeric(total.time, units = "secs")))
  return(results)
  
}      
