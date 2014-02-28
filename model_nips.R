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

#prior tau
prior_tau<-function(tau){ # tau[,i] i-th topic
  d<-sum(dgamma(tau, a0,b0,log=TRUE))
  return(d)
}

# update tau by time slice
pro_tau5<-function(N){
  cand<-matrix(data=0, nrow=slice, ncol=topic)
  cand[N,]<-rnorm(topic,0,0.1)
  #while(max(cand[N,])>700){
  #cand[N,]<-rnorm(topic,0,0.01)
  #}
  #}
  return(cand)
}
lb5<-function(beta,beta0,tau,j){# tau[j,] 
  sum<-0
  idx<-(j-1)*topic
  for(i in 1:topic){
    idx<-idx+1
    sum<-sum+sum(dnorm(beta[idx,], beta0[i,],1/tau[i],TRUE))    
  }
  return(sum)
}
ls5<-function(tau, seta,j){#tau[j,]
  #m<-vector("numeric", length(doc[j]))
  m<-vector("numeric",length(g[[j]]))
  #for(i in 1:doc[j]){ # change here
  for(i in 1:length(g[[j]])){
    m[i]<-ddirichlet(seta[g[[j]][i],], tau, log=TRUE)
    if(m[i]==Inf) print(g[[j]][i])
  }
  seta_likelihood<-sum(m)
  return(seta_likelihood)
}
up_tau5<-function(N,beta,beta0,theta,tau){
  count<-0
  exp_candpre<-tau
  #postpre<-ls(exp_candpre,theta)+lb(beta,beta0,exp_candpre/5)+prior_tau(exp_candpre)
  for(j in 1:slice){
    #j<-3
    postpre<-ls5(exp_candpre[j,],theta,j)+lb5(beta,beta0,exp_candpre[j,]/6,j)+prior_tau(exp_candpre[j,])
    for(i in 1:N){
      for(thin in 1:10){
        exp_cand<-exp(pro_tau5(j)+log(exp_candpre)) }
      post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
      a<-post-postpre
      while(a>700){
        exp_cand<-exp(pro_tau5(j)+log(exp_candpre))
        post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
        a<-post-postpre
      }
      appro<-exp(a)
      #print(appro)
      appro<-min(appro,1)
      if(runif(1)<appro){
        postpre<-post
        exp_candpre<-exp_cand
        count<-count+1
      }
      
    }
    while(count<0.17*N){
      for(i in 1:N){
        for(thin in 1:10){
          exp_cand<-exp(pro_tau5(j)+log(exp_candpre)) }
        post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
        a<-post-postpre
        while(a>700){
          exp_cand<-exp(pro_tau5(j)+log(exp_candpre))
          post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
          a<-post-postpre
        }
        appro<-exp(a)
        #print(appro)
        appro<-min(appro,1)
        if(runif(1)<appro){
          postpre<-post
          exp_candpre<-exp_cand
          count<-count+1
        }
        
      }
    }
    cat("accept for slice", j, "is ",count, "\n")
    count<-0
  }
  
  print("tau update done")
  #print(count)
  return(exp_candpre)
}

up_tau6<-function(N,beta,beta0,theta,tau){
  count<-0
  exp_candpre<-tau
  #postpre<-ls(exp_candpre,theta)+lb(beta,beta0,exp_candpre/5)+prior_tau(exp_candpre)
  for(j in 1:1){
    #j<-3
    postpre<-ls5(exp_candpre[j,],theta,j)+lb5(beta,beta0,exp_candpre[j,]/6,j)+prior_tau(exp_candpre[j,])
    for(i in 1:N){
      for(thin in 1:10){
        exp_cand<-exp(pro_tau5(j)+log(exp_candpre)) }
      post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
      a<-post-postpre
      while(a>700){
        exp_cand<-exp(pro_tau5(j)+log(exp_candpre))
        post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
        a<-post-postpre
      }
      appro<-exp(a)
      #print(appro)
      appro<-min(appro,1)
      if(runif(1)<appro){
        postpre<-post
        exp_candpre<-exp_cand
        count<-count+1
      }
      
    }
    while(count<0.17*N){
      for(i in 1:N){
        for(thin in 1:10){
          exp_cand<-exp(pro_tau5(j)+log(exp_candpre)) }
        post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
        a<-post-postpre
        while(a>700){
          exp_cand<-exp(pro_tau5(j)+log(exp_candpre))
          post<-ls5(exp_cand[j,],theta,j)+lb5(beta,beta0,exp_cand[j,]/6,j)+prior_tau(exp_cand[j,])
          a<-post-postpre
        }
        appro<-exp(a)
        #print(appro)
        appro<-min(appro,1)
        if(runif(1)<appro){
          postpre<-post
          exp_candpre<-exp_cand
          count<-count+1
        }
        
      }
    }
    cat("accept for slice", j, "is ",count, "\n")
    count<-0
  }
  
  print("tau update done")
  #print(count)
  return(exp_candpre)
}


#update Z

# get phi function
get_phi<-function(beta){
  phi<-matrix(data=NA, nrow=(topic*slice), ncol=words)
  for(i in 1:(topic*slice)){
    phi[i,]<-exp(beta[i,])/sum(exp(beta[i,]))
  }
  return(phi)
}

up_z<-function(beta,theta,dtm){
  z<-list()
  Nzv<-matrix(data=0, nrow=slice*topic,ncol=words)
  Ndz<-matrix(data=0, nrow=docs, ncol=topic)
  phi<-get_phi(beta)
  for(i in 1:slice){
    #for(j in 1:doc[i]){
    for(j in 1:length(g[[i]])){
      zrow<-which(dtm$i==g[[i]][j]) # g[[i]][j] is the doc index
      #print(zrow)
      theta_dt<-theta[g[[i]][j],] # get the topic proportation for the doc
      zterms<-dtm$j[zrow]
      topic_index<-topic*(i-1)+1 #the first topic index for i-th slice
      phi_tv<-phi[topic_index:(topic_index+topic-1), zterms] # subset of phi 
      #phi_tv<-phi[topic_index:(topic_index+7), zterms]
      zdw<-theta_dt*phi_tv
      if(length(zrow)<=1){
        zdw<-zdw/sum(zdw)
      }
      else{
        norm_zdw<-zdw/colSums(zdw)[col(zdw)] # normlization
      }
      
      #print(j)
      for(v in 1:length(zrow)){
        term_freq<-dtm$v[zrow][v]
        topic_sample<-sample(c(1:topic),term_freq, TRUE, norm_zdw[,v]) #sample topic assignment
        for(t in 1:term_freq){
          Nzv[((i-1)*topic+topic_sample[t]),zterms[v]]<-Nzv[((i-1)*topic+topic_sample[t]),zterms[v]]+1 # zterms[v], i-th word
          Ndz[g[[i]][j],topic_sample[t]]<-Ndz[g[[i]][j],topic_sample[t]]+1 # g[[i]][j]: doc index, topic_sample[t]: topic index
        }
      }
    }
  }
  z$Nzv<-Nzv
  z$Ndz<-Ndz
  print("z update done")
  return(z)
}

#update beta by time slice
pro_beta5<-function(j){
  cand<-matrix(data=0, nrow=slice*topic, ncol=words)
  for(i in 1:topic){
    idx<-(j-1)*topic+i
    cand[idx,]<-rnorm(words,0,1)
    #while(max(cand[i,])>700){
    #cand[idx,]<-rnorm(words,0,1)
    #}
  }
  return(cand)
}

lw5<-function(Nzv, phi, j){
  sum<-0
  for(i in 1:topic){
    idx<-(j-1)*topic+i # j is time slice
    sum<-sum+sum(Nzv[idx,]*log(phi[idx,]))
  }
  return(sum)
}

up_beta5<-function(N,beta,beta0,theta,tau,Nzv){ # tau is devided by 6
  count<-0
  exp_candpre<-beta
  #chain<-vector("numeric",N)
  for(j in 1:slice){
    
    phi<-get_phi(exp_candpre)
    postpre<-lb5(exp_candpre,beta0,tau[j,],j)+lw5(Nzv,phi,j)
    #postpre<-lb1(exp_candpre,beta0[j,],tau[,j],j)+lw(Nzv,phi)
    for(i in 1:N){
      exp_cand<-pro_beta5(j)+exp_candpre
      phi<-get_phi(exp_cand)
      post<-lb5(exp_cand,beta0,tau[j,],j)+lw5(Nzv,phi,j)
      #post<-lb1(exp_cand,beta0[j,],tau[,j],j)+lw(Nzv,phi)
      a<-post-postpre
      #print(str(a))
      while(a>700){
        exp_cand<-pro_beta5(j)+exp_candpre
        phi<-get_phi(exp_cand)
        post<-lb2(exp_cand,beta0,tau[j,],j)+lw5(Nzv,phi,j)
        #post<-lb1(exp_cand,beta0[j,],tau[,j],j)+lw(Nzv,phi)
        a<-post-postpre
      }
      appro<-exp(a)
      print(appro)
      appro<-min(appro,1)
      if(runif(1)<appro){
        postpre<-post
        exp_candpre<-exp_cand
        count<-count+1
      }
      #chain[i]<-exp_candpre[6,1]
    }
    cat("beta accept for topic", j ,"is", count, "\n")
    count<-0
  }
  
  print("beta update done")
  #print(count)
  return(exp_candpre)
  #return(chain)
}

#update theta

up_theta<-function(tau,Ndz){
  theta<-matrix(data=0, nrow=docs, ncol=topic)
  for(i in 1:slice){
    #for(j in 1:doc[i]){
    for(j in 1:length(g[[i]])){
      alpha<-Ndz[g[[i]][j],]+tau[i,] # g[[i]][j] the doc index
      theta[g[[i]][j],]<-rdirichlet(1, alpha)
    }
  }
  print("theta update done")
  return(theta)
}

#model

model<-function(N, dtm){
  ptm <- proc.time()
  result<-list()
  beta0<-get_beta0()
  tau<-get_tau()
  theta<-get_theta(tau)
  beta<-get_beta(beta0,tau)
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    beta0<-up_beta0(beta,tau)
    tau<-up_tau5(10,beta,beta0,theta,tau)
    z<-up_z(beta,theta,dtm)
    theta<-up_theta(tau,z$Ndz)
    phi<-get_phi(beta)
    beta<-up_beta5(10,beta,beta0,theta,tau/6,z$Nzv)
    #if(i %% 10 == 0){
    # betalist[[i]]<-beta
    #  thetalist[[i]]<-theta
    #}    
  }
  result$beta0<-beta0
  result$beta<-beta
  result$theta<-theta
  result$tau<-tau
  #result$list1<-betalist
  #result$list2<-thetalist
  print(proc.time() - ptm)
  return(result)
}

# test tau
model_tau<-function(N, dtm){
  ptm <- proc.time()
  result<-list()
  beta0<-get_beta0()
  tau<-get_tau()
  theta<-get_theta(tau)
  beta<-get_beta(beta0,tau)
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    tau<-up_tau6(10,beta,beta0,theta,tau)  
  }
  result$beta0<-beta0
  result$beta<-beta
  result$theta<-theta
  result$tau<-tau
  print(proc.time() - ptm)
  return(result)
}