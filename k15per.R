per_slda15<-vector("numeric",length(g))
for(per_i in 1:length(g)){
  per_slda15[per_i]<-perplexity(slda15[[per_i]],dtm_test[g[[per_i]],])
}