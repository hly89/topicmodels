per_slda25<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_slda25[per_i]<-perplexity(slda25[[per_i]],dtm_test[g[[per_i]],])
}