dtm<-dtm[row_sums(dtm)>0,];

#LAD
k<-30; # number of topics
seed<-2000;

gibbs<-LDA(dtm_training, k=k, method="Gibbs", control=list(seed=seed, burnin=1000, thin=100, iter=2000));