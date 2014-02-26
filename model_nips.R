# generate tau
get_tau<-function(){ 
  tau<-matrix(data=NA, nrow= slice, ncol=topic) # matrix slice*topic
  for (i in 1:slice){
    tau[i,]<-rgamma(topic, shape=a0, rate=b0)   
    while(min(tau[i,])<10^-4) {
      tau[i,]<-rgamma(topic, shape=a0, rate=b0)
    }
    
  }
  return(tau)
}

#generate beta0
get_beta0<-function(){
  beta0<-matrix(data=NA, nrow= topic, ncol=words)
  for(i in 1:topic){
    beta0[i,]<-rnorm(1:words, mean=0, sd=1)
    #beta0[i,]<-rlnorm(1:words, mean=0, sd=log(3))
  }
  return(beta0)
}

#generate theta
get_theta<-function(tau){
  theta<-matrix(data=NA, nrow=docs, ncol=topic)
  for (i in 1:slice){
    #for(j in 1:doc[i]){ 
    for(j in 1:length(g[[i]])){
      theta[g[[i]][j],]<-rdirichlet(1, tau[i,])
      while(min(theta[g[[i]][j],])<=0){
        theta[g[[i]][j],]<-rdirichlet(1, tau[i,])
      }
    }
  }
  return(theta)
}

#generate beta
get_beta<-function(beta0, tau){
  beta<-matrix(data=NA, nrow=topic*slice, ncol=words)
  for(i in 1:topic){
    for(j in 1:slice){
      idx<-(j-1)*topic+i
      #beta[idx,]<-rnorm(words, beta0[i,],1/tau[j,i])
      beta[idx,]<-rnorm(words, beta0[i,],6/tau[j,i])
    }
  }
  return(beta)
}

#update beta0
up_beta0<-function(beta, tau){
  b0<-matrix(data=NA, nrow=topic, ncol=words)
  for(i in 1:topic){
    for(j in 1:words){
      u<-0
      for(s in 1:slice){
        idx_beta<-(s-1)*topic+i
        u<-u+beta[idx_beta,j]*tau[s,i]
      }
      t<-sum(tau[,i])
      mean<-u/(1/(sd0)^2+t)
      sd<-sqrt(1/(1/(sd0^2)+t))
      b0[i,j]<-rnorm(1,mean=mean,sd=sd)
    }
  }
  print("beta0 update done")
  return(b0)
}


