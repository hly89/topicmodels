number_topic<-c(5,10,15,20,25,35,40,50,60)
gibbs_lda<-list()
for(gibbs_idx in 1:length(number_topic)){
  gibbs_lda[[gibbs_idx]]<-LDA(dtm_training, k=number_topic[gibbs_idx], method="Gibbs", control=list(seed=seed, burnin=1000, thin=100, iter=2000));
}

per_lda5<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda5[per_i]<-perplexity(gibbs_lda[[1]],dtm_test[g[[per_i]],])
}

per_lda10<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda10[per_i]<-perplexity(gibbs_lda[[2]],dtm_test[g[[per_i]],])
}

per_lda15<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda15[per_i]<-perplexity(gibbs_lda[[3]],dtm_test[g[[per_i]],])
}

per_lda20<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda20[per_i]<-perplexity(gibbs_lda[[4]],dtm_test[g[[per_i]],])
}

per_lda25<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda25[per_i]<-perplexity(gibbs_lda[[5]],dtm_test[g[[per_i]],])
}

per_lda35<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda35[per_i]<-perplexity(gibbs_lda[[6]],dtm_test[g[[per_i]],])
}

per_lda40<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda40[per_i]<-perplexity(gibbs_lda[[7]],dtm_test[g[[per_i]],])
}

per_lda50<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda50[per_i]<-perplexity(gibbs_lda[[8]],dtm_test[g[[per_i]],])
}

per_lda60<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_lda60[per_i]<-perplexity(gibbs_lda[[9]],dtm_test[g[[per_i]],])
}