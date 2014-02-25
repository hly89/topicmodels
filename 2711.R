#preprocess for time stamp
time<-apply(time,1, function(x) substr(x,1,4))
words<-words$words

#get group info
timeidx<-data.frame(table(time))$time
timeidx<-as.character(timeidx)
g<-list()
for(idx in 1:length(timeidx)){
  g[[idx]]<-which(time==timeidx[idx])
}

# get the new dtm
txt <- system.file("texts", "txt", package = "tm")
ovid <- Corpus(DirSource(txt, encoding = "UTF-8"), readerControl = list(language = "lat"))
l <- DocumentTermMatrix(ovid)
dtm<-l
dtm$i<-as.integer(i$i) 
dtm$j<-as.integer(j$j) 
dtm$v<-as.numeric(v$v)
dtm$nrow<-as.integer(length(time))
dtm$ncol<-as.integer(length(words))
dtm$dimnames$Docs<-as.character(c(1:length(time)))
dtm$dimnames$Terms<-as.character(words)

#calculate the tf-idf
tfidf<-function(dtm){
  term_tfidf<-tapply(dtm$v/row_sums(dtm)[dtm$i],dtm$j,mean)*log2(nDocs(dtm)/col_sums(dtm>0)); # row_sums from pacakage "slam"
  return(term_tfidf);
}

#filter by tf-idf
tf<-tfidf(dtm)
summary(tf);
dtm<-dtm[,tf>=0.039];
dtm<-dtm[row_sums(dtm)>0,];

# LDA global
#LAD
k<-30; # number of topics
seed<-2000;
gibbs<-LDA(dtm, k=k, method="Gibbs", control=list(seed=seed, burnin=1000, thin=100, iter=2000));
post<-posterior(gibbs);  
topics<-post$topics;# topics propotion for each doc
terms<-post$terms; # terms for 30 topics
Terms<-terms(gibbs,10);

slda<-list() # the results of the gibbs sampling 

# set the number of topics for each week is 13
for(t in 1:length(g)){
  #slda[[t]]<-LDA(dtm[g[[t]],], k=30, method="Gibbs", control=list(seed=as.integer(Sys.time()), burnin=2000, thin=100, iter=2000));
  slda[[t]]<-LDA(dtm[g[[t]],], k=30, method="Gibbs", control=list(seed=seed, burnin=2000, thin=100, iter=2000));
}

Te<-list()# to get the terms for each topic in each period
for (m in 1:N){
  Te[[m]]<-terms(g[[m]],10)
}
To<-list()# to find the top topic for each period
for (o in 1:N){
  To[[o]]<-topics(g[[o]],1)
}
order<-c(1:13)
terms_lda<-list()# get the terms distribution for each topic in each period
for(ter in 1:N){
  post_lda<-posterior(g[[ter]])
  terms_lda[[ter]]<-post_lda$terms
  #print(str(terms_lda[[ter]]))
  topicrate_lda<-data.frame(table(To[[ter]]))
  topic_order_lda<-order(-topicrate_lda$Freq)
  topic_order_lda<-c(order[topic_order_lda], order[-topic_order_lda])
  terms_lda[[ter]]<-terms_lda[[ter]][topic_order_lda,]
  Te[[ter]]<-Te[[ter]][,topic_order_lda]
}

slda[[7]]<-LDA(dtm[g[[7]],], k=30, method="Gibbs", control=list(seed=seed, burnin=2000, thin=100, iter=2000));