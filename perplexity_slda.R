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

#seperate LDA

p1<-0
p2<-sum(dtm_test)
for(i in 1:17){
  post<-posterior(slda[[i]])
  m<-length(group_test[[i]]) # doc number in slice i
  for(j in 1:m){ # for each doc
    p1<-p1+log_p(group_test[[i]][j],dtm_test,post$topics,post$terms)
  }
}

p_slda<-exp(-p1/p2)