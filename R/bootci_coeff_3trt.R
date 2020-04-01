##################################################################################
#### Bootstrap for computing the CI for coefficients alpha1, alpha2 and beta3 ####
##################################################################################

bootci_coeff_3trt <- function(NRT1, NRT2, t.seq, mediator, outcome, t.est, original.coeff, boot.sample = 1000)
{
  original.alpha1<-(original.coeff$hat.alpha1)
  original.alpha2<-(original.coeff$hat.alpha2)
  original.beta3<-(original.coeff$hat.beta3)
  
  N=length(NRT1)
  x=mediator
  y=outcome
  
  storage.boot1=matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot2=matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot3=matrix(0, nrow=boot.sample, ncol=length(t.est))
  
  start.time <- Sys.time()
  print(paste("Beginning bootstrap for coefficient CIs."))
  
  for(c1 in 1:boot.sample)
  {
    set.seed(c1^2)
    if (c1 < boot.sample) {
      cat(sprintf("Boostrapping iteration %03d", c1), " \r")
    } else {
      print(sprintf("Boostrapping iteration %03d", c1))
    }
    
    index.sample<-sample(1:N, N, replace=TRUE)
    x.boot=x[,index.sample]
    y.boot=y[,index.sample]
    NRT1.boot=NRT1[index.sample]
    NRT2.boot=NRT2[index.sample]
    deltat.boot=max(diff(t.seq))/2
    
    result.boot<-tvmcurve_3trt(NRT1.boot, NRT2.boot, t.seq, x.boot, y.boot, t.est)
    storage.boot1[c1,]<-result.boot$hat.alpha1
    storage.boot2[c1,]<-result.boot$hat.alpha2
    storage.boot3[c1,]<-result.boot$hat.beta3
  }
  
  orig.se1.all<-apply(storage.boot1, 2, sd)
  orig.se2.all<-apply(storage.boot2, 2, sd)
  orig.se3.all<-apply(storage.boot3, 2, sd)
  
  t.length<-dim(storage.boot1)[2]
  lower1<-rep(0,t.length)
  upper1<-rep(0,t.length)
  lower2<-rep(0,t.length)
  upper2<-rep(0,t.length)
  lower3<-rep(0,t.length)
  upper3<-rep(0,t.length)
  
  for(c2 in 1:(dim(storage.boot1)[2])) # run through different time points 
  {
    temp1<-storage.boot1[,c2]
    temp2<-storage.boot2[,c2]
    temp3<-storage.boot3[,c2]
    
    orig.est1<- original.alpha1[c2]
    orig.est2<- original.alpha2[c2]
    orig.est3<- original.beta3[c2]
    
    ###percentile bootstrap method 
    upper1[c2]=quantile(temp1,0.975)
    lower1[c2]=quantile(temp1,0.025)
    
    upper2[c2]=quantile(temp2,0.975)
    lower2[c2]=quantile(temp2,0.025)
    
    upper3[c2]=quantile(temp3,0.975)
    lower3[c2]=quantile(temp3,0.025)
  }
  
  coeff_CI <- list(alw1=lower1, aup1=upper1, alw2=lower2, aup2=upper2, blw3=lower3, bup3=upper3)
  
  end.time <- Sys.time()
  total.time <- end.time - start.time
  print(sprintf("Process complete. Elapsed time = %.3f secs", as.numeric(total.time, units = "secs")))
  
  return(coeff_CI)
  
}

