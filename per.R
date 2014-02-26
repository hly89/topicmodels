txt <- system.file("texts", "txt", package = "tm")
ovid <- Corpus(DirSource(txt, encoding = "UTF-8"), readerControl = list(language = "lat"))
l <- DocumentTermMatrix(ovid)

# get the new dtm
dtm$i<-as.integer(matrix$V1) # number of documents is 2340
dtm$j<-as.integer(matrix$V2) # number of terms is 21839
dtm$v<-as.numeric(matrix$V3)
dtm$nrow<-as.integer(2340)
dtm$ncol<-as.integer(21839)
dtm$dimnames$Docs<-as.character(c(1:2340))
dtm$dimnames$Terms<-as.character(words$V1)

#calculate the tf-idf
tfidf<-function(dtm){
  term_tfidf<-tapply(dtm$v/row_sums(dtm)[dtm$i],dtm$j,mean)*log2(nDocs(dtm)/col_sums(dtm>0));
  return(term_tfidf);
}

#filter by tf-idf
tf<-tfidf(dtm)
summary(tf);
dtm<-dtm[,tf>=0.039];
dtm<-dtm[row_sums(dtm)>0,];

k<-15
gibbs_20<-LDA(l, k=k, method="Gibbs", control=list(seed=as.integer(Sys.time()), burnin=1000, thin=100, iter=2000))
vem<-LDA(l, k=k, method="VEM", control=list(seed=as.integer(Sys.time())))


# compute the perplexity
perplexity_newmodel3<-vector("numeric",2000/10)
count<-1
for(i in 1:2000){
  #count<-1
  if(i %% 10 == 0){
    tzbeta<-tz$list1[[i]]
    tztheta<-tz$list2[[i]]
    tzphi<-get_phi(tzbeta)
    p1<-0
    p2<-docs*300
    for(i in 1:slice){
      m<-doc[i] # doc number in slice i
      for(j in 1:m){ # for each doc
        p1<-p1+log_p(g[[i]][j],dtm1,tztheta,tzphi[((i-1)*15+1):(i*15),])
      }
    }
    #print(exp(-p1/p2))
    perplexity_newmodel3[count] = exp(-p1/p2)
    count<-count+1
  }
}

perplexity_newmodel10<-vector("numeric",11000/10)
count<-1
for(i in 1:11000){
  #count<-1
  if(i %% 10 == 0){
    tzbeta<-tz10$list1[[i]]
    tztheta<-tz10$list2[[i]]
    tzphi<-get_phi(tzbeta)
    p1<-0
    p2<-docs*300
    for(i in 1:slice){
      m<-doc[i] # doc number in slice i
      for(j in 1:m){ # for each doc
        p1<-p1+log_p(g[[i]][j],dtm1,tztheta,tzphi[((i-1)*15+1):(i*15),])
      }
    }
    #print(exp(-p1/p2))
    perplexity_newmodel10[count] = exp(-p1/p2)
    count<-count+1
  }
}

log_p<-function(doc_no,dtm,theta,phi){ # dtm is test dataset, phi in one time slice
  idx<-which(dtm$i == doc_no)
  sum<-0
  for(i in 1:length(idx)){
    n_dt<-dtm$v[idx[i]] # number of times for the term t in doc d
    t_d<-dtm$j[idx[i]] # term t index
    sum<-sum+n_dt*log(sum(phi[,t_d]*theta[doc_no,]))
  }
  return(sum)
}

log_p1<-function(doc_no,dtm,theta,phi){ # dtm is test dataset, phi in one time slice
  idx<-which(dtm$i == doc_no)
  if(doc_no%%20 == 0){
    theta_idx<-20
  }
  else{
    theta_idx<-doc_no%%20
  }
  sum<-0
  for(i in 1:length(idx)){
    n_dt<-dtm$v[idx[i]] # number of times for the term t in doc d
    t_d<-dtm$j[idx[i]] # term t index
    sum<-sum+n_dt*log(sum(phi[,t_d]*theta[theta_idx,]))
  }
  return(sum)
}

p1<-0
p2<-docs*300
for(i in 1:slice){
  m<-doc[i] # doc number in slice i
  for(j in 1:m){ # for each doc
    p1<-p1+log_p(g[[i]][j],dtm1,tz$theta,tzphi[((i-1)*15+1):(i*15),])
  }
}

p<-exp(-p1/p2)

# perplexity for all docs in LDA(gibbs)
p1<-0
p2<-docs*300
for(i in 1:slice){
  m<-doc[i] # doc number in slice i
  for(j in 1:m){ # for each doc
    p1<-p1+log_p(g[[i]][j],dtm1,post$topics,post$terms)
  }
}

p<-exp(-p1/p2)

# vem
post_vem<-posterior(vem)
p1<-0
p2<-docs*300
for(i in 1:slice){
  m<-doc[i] # doc number in slice i
  for(j in 1:m){ # for each doc
    p1<-p1+log_p(g[[i]][j],dtm1,post_vem$topics,post_vem$terms)
  }
}

p<-exp(-p1/p2) # 293

perplexity_newmodel5<-vector("numeric",3500/10)
count<-1
for(i in 1:3500){
  #count<-1
  if(i %% 10 == 0){
    tzbeta<-tz$list1[[i]]
    tztheta<-tz$list2[[i]]
    tzphi<-get_phi(tzbeta)
    p1<-0
    p2<-docs*300
    for(i in 1:slice){
      m<-doc[i] # doc number in slice i
      for(j in 1:m){ # for each doc
        p1<-p1+log_p(g[[i]][j],dtm1,tztheta,tzphi[((i-1)*15+1):(i*15),])
      }
    }
    #print(exp(-p1/p2))
    perplexity_newmodel5[count] = exp(-p1/p2)
    count<-count+1
  }
}#min:283,1566

perplexity_newmodel6<-vector("numeric",6000/10)
count<-1
for(i in 1:6000){
  #count<-1
  if(i %% 10 == 0){
    tzbeta<-tz6$list1[[i]]
    tztheta<-tz6$list2[[i]]
    tzphi<-get_phi(tzbeta)
    p1<-0
    p2<-docs*300
    for(i in 1:slice){
      m<-doc[i] # doc number in slice i
      for(j in 1:m){ # for each doc
        p1<-p1+log_p(g[[i]][j],dtm1,tztheta,tzphi[((i-1)*15+1):(i*15),])
      }
    }
    #print(exp(-p1/p2))
    perplexity_newmodel6[count] = exp(-p1/p2)
    count<-count+1
  }
}#min:

#seperate LDA
for(t in 1:5){
  lda_sep[[t]]<-LDA(l[g[[t]],], k=15, method="Gibbs", control=list(seed=as.integer(Sys.time()), burnin=2000, thin=100, iter=2000));
}
p1<-0
p2<-docs*300
for(i in 1:slice){
  post<-posterior(lda_sep[[i]])
  m<-doc[i] # doc number in slice i
  for(j in 1:m){ # for each doc
    p1<-p1+log_p(g[[i]][j],dtm1,post$topics,post$terms)
  }
}

p<-exp(-p1/p2)

#perplexity for slice 1(first 20 docs)
p1<-0
p2<-20*300
for(i in 1:20){
  p1<-p1+log_p(g[[1]][i],dtm1,post$topics,post$terms)
}
p<-exp(-p1/p2) # p:275,8262

perplexity_slice1<-vector("numeric",2000/10)
count<-1
for(i in 1:2000){
  #count<-1
  if(i %% 10 == 0){
    tzbeta<-tz$list1[[i]]
    tztheta<-tz$list2[[i]]
    tzphi<-get_phi(tzbeta)
    p1<-0
    p2<-20*300
    #for(i in 1:slice){
      #m<-doc[i] # doc number in slice i
      for(j in 1:20){ # for each doc
        p1<-p1+log_p(g[[1]][j],dtm1,tztheta,tzphi[1:15,])
      }
    #}
    #print(exp(-p1/p2))
    perplexity_slice1[count] = exp(-p1/p2)
    count<-count+1
  }
}# min:370,5495

post1<-posterior(lda_sep[[1]])
p1<-0
p2<-20*300
for(i in 1:20){
  p1<-p1+log_p(g[[1]][i],dtm1,post1$topics,post1$terms)
}
p<-exp(-p1/p2) #p:291,5035


#perplexity for slice 2(20:40 docs)
p1<-0
p2<-20*300
for(i in 1:20){
  p1<-p1+log_p(g[[2]][i],dtm1,post$topics,post$terms)
}

p<-exp(-p1/p2)# p:195,8815

perplexity_slice2<-vector("numeric",2000/10)
count<-1
for(i in 1:2000){
  #count<-1
  if(i %% 10 == 0){
    tzbeta<-tz$list1[[i]]
    tztheta<-tz$list2[[i]]
    tzphi<-get_phi(tzbeta)
    p1<-0
    p2<-20*300
    #for(i in 1:slice){
    #m<-doc[i] # doc number in slice i
    for(j in 1:20){ # for each doc
      p1<-p1+log_p(g[[2]][j],dtm1,tztheta,tzphi[16:30,])
    }
    #}
    #print(exp(-p1/p2))
    perplexity_slice2[count] = exp(-p1/p2)
    count<-count+1
  }
}# min:261,7484

post2<-posterior(lda_sep[[2]])
p1<-0
p2<-20*300
for(i in 1:20){
  p1<-p1+log_p1(g[[2]][i],dtm1,post2$topics,post2$terms)
  print(p1)
}
p<-exp(-p1/p2)#p:156,9457