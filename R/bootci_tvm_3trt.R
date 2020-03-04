## function to calculate the time-varying-mediation-curve and standard deviation
bootci_tvm_3trt<-function(boot.sample, orig.data, t.est)
{
  N=length(orig.data$NRT1)
  x=orig.data$x
  y=orig.data$y
  NRT1=orig.data$NRT1
  NRT2=orig.data$NRT2
  t.seq<-orig.data$t.seq
  deltat=max(diff(t.seq))/2
  result.org <- tvmcurve_3trt(NRT1, NRT2, t.seq, x, y, t.est)
  # result.org<-tvmcurve_3trt_sd(NRT1, NRT2, t.seq, x, y, t.est) 
  original.mediation1<-(result.org$hat.mediation1)
  original.mediation2<-(result.org$hat.mediation2)
  # orig.sd1.all<-result.org$sd1
  # orig.sd2.all<-result.org$sd2
  
  # orig.se<-tvmcurve_boot_se(orig.data, t.est, boot.sample)
  # orig.se1.all<-orig.se$boot.se1
  # orig.se2.all<-orig.se$boot.se2
  
  storage.boot1=matrix(0, nrow=boot.sample, ncol=length(t.est))
  storage.boot2=matrix(0, nrow=boot.sample, ncol=length(t.est))
  
  print(paste("Bootstrap starting for confidence interval computation."))
  
  for(k1 in 1:boot.sample)
  {
    set.seed(k1^2)
    print(k1)
    index.sample<-sample(1:N, N, replace=TRUE)
    x.boot=x[,index.sample]
    y.boot=y[,index.sample]
    NRT1.boot=NRT1[index.sample]
    NRT2.boot=NRT2[index.sample]
    deltat=max(diff(t.seq))/2
    
    result.boot<-tvmcurve_3trt(NRT1.boot, NRT2.boot, t.seq, x.boot, y.boot, t.est) 
    storage.boot1[k1,]<-result.boot$hat.mediation1
    storage.boot2[k1,]<-result.boot$hat.mediation2
  }
  
  orig.se1.all<-apply(storage.boot1, 2, sd)
  orig.se2.all<-apply(storage.boot2, 2, sd)
  # result.sd1<-apply(storage.boot1, 2, sd)
  # result.sd2<-apply(storage.boot2, 2, sd)
  # 
  t.length<-dim(storage.boot1)[2]
  lower1<-rep(0,t.length)
  upper1<-rep(0,t.length)
  lower2<-rep(0,t.length)
  upper2<-rep(0,t.length)
  
  for(k2 in 1:(dim(storage.boot1)[2])) # run through different time points 
  {
    temp1<-storage.boot1[,k2]
    temp2<-storage.boot2[,k2]
    
    orig.est1<- original.mediation1[k2]
    orig.est2<- original.mediation2[k2]
    
    ###percentile bootstrap method 
    upper1[k2]=quantile(temp1,0.975)
    lower1[k2]=quantile(temp1,0.025)
    
    upper2[k2]=quantile(temp2,0.975)
    lower2[k2]=quantile(temp2,0.025)
  }
  
  ## The line below is commented till orig.sd1.all and orig.sd2.all from tvmcurve_3trt_sd.R can be computed
  # list(plw1=lower1, pup1=upper1, plw2=lower2, pup2=upper2,  orig.sd1.all=orig.sd1.all, orig.sd2.all=orig.sd2.all, orig.se1.all=orig.se1.all, orig.se2.all=orig.se2.all, orig.mediation1=original.mediation1, orig.mediation2=original.mediation2) 
  
  list(plw1=lower1, pup1=upper1, plw2=lower2, pup2=upper2, orig.se1.all=orig.se1.all, orig.se2.all=orig.se2.all, orig.mediation1=original.mediation1, orig.mediation2=original.mediation2) 
  
}      
