dtm<-dtm[row_sums(dtm)>0,];
dtm_training<-dtm_training[row_sums(dtm_training)>0,];
dtm_test<-dtm_test[row_sums(dtm_test)>0,];
#LAD
k<-30; # number of topics
seed<-2000;

gibbs<-LDA(dtm_training, k=k, method="Gibbs", control=list(seed=seed, burnin=1000, thin=100, iter=2000));


slda<-list() # the results of the gibbs sampling 

# set the number of topics for each week is 13
for(t in 1:length(g)){
  #slda[[t]]<-LDA(dtm[g[[t]],], k=30, method="Gibbs", control=list(seed=as.integer(Sys.time()), burnin=2000, thin=100, iter=2000));
  slda[[t]]<-LDA(dtm_training[g[[t]],], k=30, method="Gibbs", control=list(seed=seed, burnin=2000, thin=100, iter=2000));
}

per_slda<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_slda[per_i]<-perplexity(slda[[per_i]],dtm_test[g[[per_i]],])
}

per_lda<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda[per_i]<-perplexity(gibbs,dtm_test[g[[per_i]],])
}