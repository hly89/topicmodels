ptm <- proc.time()
slda50<-list() # the results of the gibbs sampling 

# set the number of topics for each week is 5
for(t in 1:length(g)){
  #slda[[t]]<-LDA(dtm[g[[t]],], k=30, method="Gibbs", control=list(seed=as.integer(Sys.time()), burnin=2000, thin=100, iter=2000));
  cat("Iteration: ", t,"\n")
  slda50[[t]]<-LDA(dtm_training[g[[t]],], k=50, method="Gibbs", control=list(seed=seed, burnin=2000, thin=100, iter=2000));
}
print(proc.time() - ptm)