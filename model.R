# generate tau
get_tau<-function(){
  tau<-matrix(data=NA, nrow= slice, ncol=topic)
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
    for(j in 1:doc[i]){
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


#update tau

#proposal function

pro_tau1<-function(N){
  cand<-matrix(data=0, nrow=slice, ncol=topic)
  #for(i in 1:2){
    cand[,N]<-rnorm(slice,0,0.1)
    while(max(cand[,N])>700){
      cand[,N]<-rnorm(slice,0,0.1)
    }
  #}
  return(cand)
}

#likelihood
ls<-function(tau, seta){
  t<-vector("numeric", slice)
  for(i in 1:slice){
    m<-vector("numeric", doc[i])
    for(j in 1:doc[i]){
      m[j]<-ddirichlet(seta[g[[i]][j],], tau[i,], log=TRUE)
      if(m[j]==Inf) print(g[[i]][j])
    }
    t[i]<-sum(m)
  }
  seta_likelihood<-sum(t)
  return(seta_likelihood)
}

lb<-function(beta,beta0,tau){
  idx<-1
  sum<-0
  for(i in 1:slice){
    for(j in 1:topic){
      for(v in 1:words){
        sum<-sum+dnorm(beta[idx,v], beta0[j,v],1/tau[i,j],TRUE)
      }
      idx<-idx+1
    }
  }
  return(sum)
}

lb1<-function(beta,beta0,tau,j){# tau[,j] j-th topic beta0[j,]  
  sum<-0
  for(i in 1:slice){
    idx<-(i-1)*15+j
    for(v in 1:words){      
      sum<-sum+dnorm(beta[idx,v], beta0[v],1/tau[i],TRUE)
      #print(sum)
    }      
  }
  return(sum)
}

lb2<-function(beta,beta0,tau,j){# tau[,j] j-th topic beta0[j,]  
  sum<-0
  for(i in 1:slice){
    idx<-(i-1)*topic+j
    sum<-sum+sum(dnorm(beta[idx,], beta0,1/tau[i],TRUE))    
  }
  return(sum)
}

#prior tau
prior_tau<-function(tau){ # tau[,i] i-th topic
  d<-sum(dgamma(tau, a0,b0,log=TRUE))
  return(d)
}


up_tau<-function(N,beta,beta0,theta,tau){
  count<-0
  exp_candpre<-tau
  postpre<-ls(exp_candpre,theta)+lb(beta,beta0,exp_candpre)+prior_tau(exp_candpre)
  for(i in 1:N){
    exp_cand<-exp(pro_tau1()+log(exp_candpre))
    post<-ls(exp_cand,theta)+lb(beta,beta0,exp_cand)+prior_tau(exp_cand)
    a<-post-postpre
    while(a>700){
      exp_cand<-exp(pro_tau1()+log(exp_candpre))
      post<-ls(exp_cand,theta)+lb(beta,beta0,exp_cand)+prior_tau(exp_cand)
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
  }
  print("tau update done")
  print(count)
  return(exp_candpre)
}

up_tau1<-function(N,beta,beta0,theta,tau){
  count<-0
  exp_candpre<-tau
  postpre<-ls(exp_candpre,theta)+lb(beta,beta0,exp_candpre)+prior_tau(exp_candpre)
  for(j in 1:topic){
    for(i in 1:N){
      exp_cand<-exp(pro_tau1(j)+log(exp_candpre))
      post<-ls(exp_cand,theta)+lb(beta,beta0,exp_cand)+prior_tau(exp_cand)
      a<-post-postpre
      while(a>700){
        exp_cand<-exp(pro_tau1(j)+log(exp_candpre))
        post<-ls(exp_cand,theta)+lb(beta,beta0,exp_cand)+prior_tau(exp_cand)
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
      
    }
    cat("accept for topic", j, "is ",count, "\n")
    count<-0
  }
  
  print("tau update done")
  #print(count)
  return(exp_candpre)
}

# change tau for lb (update tau by topic)
up_tau2<-function(N,beta,beta0,theta,tau){
  count<-0
  exp_candpre<-tau
  #postpre<-ls(exp_candpre,theta)+lb(beta,beta0,exp_candpre/5)+prior_tau(exp_candpre)
  for(j in 1:topic){
    postpre<-ls(exp_candpre,theta)+lb2(beta,beta0[j,],exp_candpre[,j]/5,j)+prior_tau(exp_candpre[,j])
    for(i in 1:N){
      exp_cand<-exp(pro_tau1(j)+log(exp_candpre))
      post<-ls(exp_cand,theta)+lb2(beta,beta0[j,],exp_cand[,j]/5,j)+prior_tau(exp_cand[,j])
      a<-post-postpre
      while(a>700){
        exp_cand<-exp(pro_tau1(j)+log(exp_candpre))
        post<-ls(exp_cand,theta)+lb2(beta,beta0[j,],exp_cand[,j]/5,j)+prior_tau(exp_cand[,j])
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
      
    }
    cat("accept for topic", j, "is ",count, "\n")
    count<-0
  }
  
  print("tau update done")
  #print(count)
  return(exp_candpre)
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
  m<-vector("numeric", length(doc[j]))

    for(i in 1:doc[j]){
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

up_tau3<-function(N,beta,beta0,theta,tau){
  count<-0
  exp_candpre<-tau
  chain<-vector("numeric",N)
  #postpre<-ls(exp_candpre,theta)+lb(beta,beta0,exp_candpre/5)+prior_tau(exp_candpre)
  #for(j in 1:topic){
  j<-2
    postpre<-ls(exp_candpre,theta)+lb2(beta,beta0[j,],exp_candpre[,j]/5,j)+prior_tau(exp_candpre[,j])
    for(i in 1:N){
      exp_cand<-exp(pro_tau1(j)+log(exp_candpre))
      post<-ls(exp_cand,theta)+lb2(beta,beta0[j,],exp_cand[,j]/5,j)+prior_tau(exp_cand[,j])
      a<-post-postpre
      while(a>700){
        exp_cand<-exp(pro_tau1(j)+log(exp_candpre))
        post<-ls(exp_cand,theta)+lb2(beta,beta0[j,],exp_cand[,j]/5,j)+prior_tau(exp_cand[,j])
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
      chain[i]<-exp_candpre[1,2]
    }
    #cat("accept for topic", j, "is ",count, "\n")
    #count<-0
  #}
  
  print("tau update done")
  print(count)
  #return(exp_candpre)
  return(chain)
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
    for(j in 1:doc[i]){
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

# update beta

#proposal for beta
pro_beta<-function(){
  cand<-matrix(data=NA, nrow=slice*topic, ncol=words)
  for(i in 1:(slice*topic)){
    cand[i,]<-rnorm(words,0,0.01)
    while(max(cand[i,])>700){
      cand[i,]<-rnorm(words,0,0.01)
    }
  }
  return(cand)
}

pro_beta1<-function(j){
  cand<-matrix(data=0, nrow=slice*topic, ncol=words)
  for(i in 1:slice){
    idx<-(i-1)*15+j
    cand[idx,]<-rnorm(words,0,0.08)
    while(max(cand[i,])>700){
      cand[idx,]<-rnorm(words,0,0.08)
    }
  }
  return(cand)
}

#likelihood 
lw<-function(Nzv, phi){
  sum<-0
  topic_index<-1
  for(i in 1:slice){
    for(j in 1:topic){
      for(v in 1:words){
        sum<-sum+Nzv[topic_index,v]*log(phi[topic_index,v])
      }
      topic_index<-topic_index+1
      #print(topic_index)
    }
  }
  return(sum)
}


# update for each topic
lw1<-function(Nzv, phi, j){
  sum<-0
  for(i in 1:slice){
    idx<-(i-1)*15+j
    sum<-sum+sum(Nzv[idx,]*log(phi[idx,]))
  }
  return(sum)
}

up_beta<-function(N,beta,beta0,tau,Nzv){
  count<-0
  exp_candpre<-beta
  for(j in 1:topic){
    phi<-get_phi(exp_candpre)
    postpre<-lb2(exp_candpre,beta0[j,],tau[,j],j)+lw1(Nzv,phi,j)
    #postpre<-lb1(exp_candpre,beta0[j,],tau[,j],j)+lw(Nzv,phi)
    for(i in 1:N){
      exp_cand<-pro_beta1(j)+exp_candpre
      phi<-get_phi(exp_cand)
      post<-lb2(exp_cand,beta0[j,],tau[,j],j)+lw1(Nzv,phi,j)
      #post<-lb1(exp_cand,beta0[j,],tau[,j],j)+lw(Nzv,phi)
      a<-post-postpre
      while(a>700){
        exp_cand<-pro_beta1(j)+exp_candpre
        phi<-get_phi(exp_cand)
        post<-lb2(exp_cand,beta0[j,],tau[,j],j)+lw1(Nzv,phi,j)
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
    }
    cat("beta accept for topic", j ,"is", count, "\n")
    count<-0
  }
  
  print("beta update done")
  #print(count)
  return(exp_candpre)
}

up_beta2<-function(N,beta,beta0,theta,tau,Nzv){
  count<-0
  exp_candpre<-beta
  chain<-vector("numeric",N)
  #for(j in 1:2){
    j<-6
    phi<-get_phi(exp_candpre)
    postpre<-lb2(exp_candpre,beta0[j,],tau[,j],j)+lw1(Nzv,phi,j)
    #postpre<-lb1(exp_candpre,beta0[j,],tau[,j],j)+lw(Nzv,phi)
    for(i in 1:N){
      exp_cand<-pro_beta1(j)+exp_candpre
      phi<-get_phi(exp_cand)
      post<-lb2(exp_cand,beta0[j,],tau[,j],j)+lw1(Nzv,phi,j)
      #post<-lb1(exp_cand,beta0[j,],tau[,j],j)+lw(Nzv,phi)
      a<-post-postpre
      while(a>700){
        exp_cand<-pro_beta1(j)+exp_candpre
        phi<-get_phi(exp_cand)
        post<-lb2(exp_cand,beta0[j,],tau[,j],j)+lw1(Nzv,phi,j)
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
      chain[i]<-exp_candpre[6,1]
    }
    cat("beta accept for topic", j ,"is", count, "\n")
    #count<-0
  #}
  
  print("beta update done")
  #print(count)
  #return(exp_candpre)
  return(chain)
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

up_beta5<-function(N,beta,beta0,theta,tau,Nzv){ # tau is devided by 8
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
    for(j in 1:doc[i]){
      alpha<-Ndz[g[[i]][j],]+tau[i,] # g[[i]][j] the doc index
      theta[g[[i]][j],]<-rdirichlet(1, alpha)
    }
  }
  print("theta update done")
  return(theta)
}

# updata beta one by one
up_betaS<-function(N,beta,beta0,tau,Nzv, phi){ 
  #ptm <- proc.time() 
  #chain<-vector("numeric", N)
  exp_candpre<-beta
  postpre<-dnorm(exp_candpre,beta0,1/tau,TRUE)+Nzv*log(phi)
  for(i in 1:N){
    exp_cand<-rnorm(1,0,0.08)+exp_candpre
    post<-dnorm(exp_cand,beta0,1/tau,TRUE)+Nzv*log(phi)
    a<-post-postpre
    while(a>700){
      exp_cand<-rnorm(1,0,0.08)+exp_candpre
      post<-dnorm(exp_cand,beta0,1/tau,TRUE)+Nzv*log(phi)
      a<-post-postpre
    }
    appro<-exp(a)
    #print(appro)
    appro<-min(appro,1)
    if(runif(1)<appro){
      postpre<-post
      exp_candpre<-exp_cand
    }
    #chain[i]<-exp_candpre
  }
  #cat("beta accept is", count, "\n")
  
  #print("beta update done")
  #print(proc.time() - ptm)
  #print(count)
  return(exp_candpre)
}

#model

model<-function(N, dtm){
  ptm <- proc.time()
  result<-list()
  #betalist<-list()
  #thetalist<-list()
  beta0<-get_beta0()
  tau<-get_tau()
  theta<-get_theta(tau)
  beta<-get_beta(beta0,tau)
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    beta0<-up_beta0(beta,tau)
    tau<-up_tau5(100,beta,beta0,theta,tau)
    z<-up_z(beta,theta,dtm)
    theta<-up_theta(tau,z$Ndz)
    phi<-get_phi(beta)
    for(t in 1:slice){
      for(j in 1:topic){
        topic_idx<-(t-1)*topic+j
        #print(topic_idx)
        for(v in 1:words){
          beta[topic_idx,v]<-up_betaS(10,beta[topic_idx,v],beta0[j,v],tau[t,j]/6,z$Nzv[topic_idx,v],phi[topic_idx,v])
          #print(v)
        }
      }
    }
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

model_beta0<-function(N, dtm,beta0){
  ptm <- proc.time()
  result<-list()
  #betalist<-list()
  #thetalist<-list()
  #beta0<-get_beta0()
  tau<-get_tau()
  theta<-get_theta(tau)
  beta<-get_beta(beta0,tau)
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    #beta0<-up_beta0(beta,tau)
    tau<-up_tau5(10,beta,beta0,theta,tau)
    z<-up_z(beta,theta,dtm)
    theta<-up_theta(tau,z$Ndz)
    phi<-get_phi(beta)
    for(t in 1:slice){
      for(j in 1:topic){
        topic_idx<-(t-1)*topic+j
        #print(topic_idx)
        for(v in 1:words){
          beta[topic_idx,v]<-up_betaS(10,beta[topic_idx,v],beta0[j,v],tau[t,j]/6,z$Nzv[topic_idx,v],phi[topic_idx,v])
          #print(v)
        }
      }
    }
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

model_1<-function(N, dtm,beta0,theta,Nzv){
  ptm <- proc.time()
  result<-list()
  thetalist<-list()
  tau<-get_tau()
  beta<-get_beta(beta0,tau)
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    tau<-up_tau5(10,beta,beta0,theta,tau)
    phi<-get_phi(beta)
    for(t in 1:slice){
      for(j in 1:topic){
        topic_idx<-(t-1)*topic+j
        #print(topic_idx)
        for(v in 1:words){
          beta[topic_idx,v]<-up_betaS(100,beta[topic_idx,v],beta0[j,v],tau[t,j]/6,Nzv[topic_idx,v],phi[topic_idx,v])
          #print(v)
        }
      }
    }
    if(i %% 10 == 0){
      thetalist[[i]]<-theta
    }
    
  }
  result$beta<-beta
  result$tau<-tau
  result$list2<-thetalist
  print(proc.time() - ptm)
  return(result)
}

model_2<-function(N, dtm,tau,theta,Nzv){
  ptm <- proc.time()
  result<-list()
  thetalist<-list()
  beta0<-get_beta0()
  beta<-get_beta(beta0,tau)
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    beta0<-up_beta0(beta,tau)
    phi<-get_phi(beta)
    for(t in 1:slice){
      for(j in 1:topic){
        topic_idx<-(t-1)*topic+j
        #print(topic_idx)
        for(v in 1:words){
          beta[topic_idx,v]<-up_betaS(100,beta[topic_idx,v],beta0[j,v],tau[t,j]/6,Nzv[topic_idx,v],phi[topic_idx,v])
          #print(v)
        }
      }
    }
    if(i %% 10 == 0){
      thetalist[[i]]<-theta
    }
    
  }
  result$beta<-beta
  result$beta0<-beta0
  result$list2<-thetalist
  print(proc.time() - ptm)
  return(result)
}

model_3<-function(N, dtm,beta0,beta,Ndz){
  ptm <- proc.time()
  result<-list()
  thetalist<-list()
  tau<-get_tau()
  theta<-get_theta(tau)
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    tau<-up_tau5(10,beta,beta0,theta,tau)
    theta<-up_theta(tau,Ndz)
    
  }
  result$tau<-tau
  result$theta<-theta
  print(proc.time() - ptm)
  return(result)
}

model_4<-function(N, dtm, beta,theta, Nzv){
  ptm <- proc.time()
  result<-list()
  beta0<-get_beta0()
  tau<-get_tau()
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    beta0<-up_beta0(beta,tau)
    tau<-up_tau5(10,beta,beta0,theta,tau)
  }
  result$beta0<-beta0
  result$tau<-tau
  print(proc.time() - ptm)
  return(result)
}

model1<-function(N, dtm, beta){
  result<-list()
  beta0<-get_beta0()
  tau<-get_tau()
  theta<-get_theta(tau)
  beta1<-get_beta(beta0,tau)
  beta1[1:3,]<-beta[1:3,]
  for(i in 1:N){
    cat("Iteration: ", i,"\n")
    beta0<-up_beta0(beta1,tau)
    tau<-up_tau5(10,beta1,beta0,theta,tau)
    z<-up_z(beta1,theta,dtm)
    theta<-up_theta(tau,z$Ndz)
    beta1<-up_beta5(10,beta1,beta0,theta,tau/5,z$Nzv)
  }
  result$beta0<-beta0
  result$beta<-beta1
  result$theta<-theta
  result$tau<-tau
  return(result)
}

#test beta0 part
modelb0<-function(N, dtm, tau, beta, theta, beta0_init){
  result<-list()
  beta0_init<-beta0_init
  beta0_pre<-beta0_init
  beta0<-get_beta0()
  disinit<-vector("numeric",N)
  dispre<-vector("numeric",N)
  for(i in 1:N){
    print("Iteration: ", i)
    beta0<-up_beta0(beta,tau)
    disinit[i]<-sqrt(sum((beta0-beta0_init)^2))
    dispre[i]<-sqrt(sum((beta0-beta0_pre)^2))
    beta0_pre<-beta0
    tau<-up_tau(10,beta,beta0,theta,tau)
    z<-up_z(beta,theta,dtm)
    theta<-up_theta(tau,z$Ndz)
    beta<-up_beta(10,beta,beta0,theta,tau,z$Nzv)
  }
  result$beta0<-beta0
  result$beta<-beta
  result$theta<-theta
  result$tau<-tau
  result$disinit<-disinit
  result$dispre<-dispre
  return(result)
}

#test tau part
modeltau<-function(N, dtm,beta,tau_init, theta){
  result<-list()
  tau_init<-tau_init
  tau_pre<-tau_init
  tau<-get_tau()
  disinit<-vector("numeric",N)
  dispre<-vector("numeric",N)
  for(i in 1:N){
    cat("Iteration: ", i, "\n")
    beta0<-up_beta0(beta,tau)    
    tau<-up_tau2(10,beta,beta0,theta,tau)
    disinit[i]<-sqrt(sum((tau-tau_init)^2))
    dispre[i]<-sqrt(sum((tau-tau_pre)^2))
    tau_pre<-tau
    z<-up_z(beta,theta,dtm)
    theta<-up_theta(tau,z$Ndz)
    beta<-up_beta(10,beta,beta0,theta,tau/5,z$Nzv)
  }
  result$beta0<-beta0
  result$beta<-beta
  result$theta<-theta
  result$tau<-tau
  result$disinit<-disinit
  result$dispre<-dispre
  return(result)
}

#test theta part
modeltheta<-function(){
  
}
