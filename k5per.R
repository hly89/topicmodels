per_slda5<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_slda5[per_i]<-perplexity(slda5[[per_i]],dtm_test[g[[per_i]],])
}