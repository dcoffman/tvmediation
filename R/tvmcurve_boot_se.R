tvmcurve_boot_se<-function(data, t.est, num.boot.sample=1000)
{
  deltat=max(diff(t.seq))/2
  storage.boot1.sub=matrix(0, nrow=num.boot.sample, ncol=length(t.est))
  storage.boot2.sub=matrix(0, nrow=num.boot.sample, ncol=length(t.est))
  N=length(data$NRT1) #number of cases
  t.seq=data$t.seq
  print("Bootstrap starting for standard error computation.")
  
  for(l in 1:num.boot.sample)
  {  
    set.seed(l^2)
    index.sample<-sample(1:N, N, replace=TRUE)
    x.boot.sub=data$x[,index.sample]
    y.boot.sub=data$y[,index.sample]
    NRT1.boot.sub=data$NRT1[index.sample]
    NRT2.boot.sub=data$NRT2[index.sample]
    result.boot<-tvmcurve_3trt(NRT1.boot.sub, NRT2.boot.sub, t.seq, x.boot.sub, y.boot.sub, t.est)
    
    storage.boot1.sub[l,]<-(result.boot$hat.mediation1)
    storage.boot2.sub[l,]<-(result.boot$hat.mediation2)
    
    print(paste("Bootstrap sample",l,"complete."))
  }
  result.sd1<-apply(storage.boot1.sub, 2, sd)
  result.sd2<-apply(storage.boot2.sub, 2, sd)
  
  list(boot.se1=result.sd1, boot.se2=result.sd2)
}
