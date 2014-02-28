number_topic<-c(5,10,15,20,25,35,40,50,60)
gibbs_lda<-list()
for(gibbs_idx in 1:length(number_topic)){
  gibbs_lda[[gibbs_idx]]<-LDA(dtm_training, k=number_topic[gibbs_idx], method="Gibbs", control=list(seed=seed, burnin=1000, thin=100, iter=2000));
}
