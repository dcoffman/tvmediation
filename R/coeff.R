###### It assumes the following quantities are known #######
#### NRT: intervention 
#### t-seq: observed t-seq 
#### x, y, deltat
coeff<-function(j, NRT1, NRT2, x, y)
{
  X.new.j<-cbind(c(NRT1, rep(0, length(NRT1))),c(NRT2, rep(0, length(NRT2))),c(rep(0, length(NRT1)), NRT1),c(rep(0, length(NRT2)), NRT2),c(rep(0, length(NRT1)), as.numeric(x[j-1,])))
  Y.new.j<-c(x[j,], y[j-1,])
  
  X.new.j=scale(X.new.j,center=TRUE, scale=FALSE)
  Y.new.j=scale(Y.new.j,center=TRUE, scale=FALSE)
  
  nomissing.x<-complete.cases(X.new.j)
  nomissing.y<-complete.cases(Y.new.j)
  
  nomissing.index <- nomissing.x * nomissing.y
  
  X.new.j=X.new.j[which(nomissing.index==1),]
  Y.new.j=Y.new.j[which(nomissing.index==1)]
  
  coeff.est<-solve(t(X.new.j)%*%(X.new.j))%*%t(X.new.j)%*%(Y.new.j)
  list(coeff.est=coeff.est, nomissing.index=nomissing.index)
}
