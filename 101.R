training_i<-c()
training_j<-c()
training_v<-c()
test_i<-c()
test_j<-c()
test_v<-c()
for(doc in 1801:1850){ # change here
  wordcounts<-dtm$v[which(dtm$i==doc)]
  termsindex<-dtm$j[which(dtm$i==doc)]
  for(word in 1: length(wordcounts)){
    if(wordcounts[word]%%2==0){
      training_i<-cbind(training_i, doc)
      test_i<-cbind(test_i, doc)
      training_j<-cbind(training_j, termsindex[word])
      test_j<-cbind(test_j, termsindex[word])
      training_v<-cbind(training_v, wordcounts[word]/2)
      test_v<-cbind(test_v, wordcounts[word]/2)
    } else if(wordcounts[word]==1&&word%%2==0){
      test_i<-cbind(test_i, doc)
      test_j<-cbind(test_j, termsindex[word])
      test_v<-cbind(test_v, 1)
    }else if(wordcounts[word]==1&&word%%2!=0){
      training_i<-cbind(training_i, doc)
      training_j<-cbind(training_j, termsindex[word])
      training_v<-cbind(training_v, 1)
    }else{
      training_i<-cbind(training_i, doc)
      training_j<-cbind(training_j, termsindex[word])
      training_v<-cbind(training_v, round(wordcounts[word]/2))
      test_i<-cbind(test_i, doc)
      test_j<-cbind(test_j, termsindex[word])
      test_v<-cbind(test_v, wordcounts[word]-round(wordcounts[word]/2))
    }
  }
}

# from 101 to 500
training_i<-as.integer(training_i)
dtm_training$i<-append(dtm_training$i,training_i)
training_j<-as.integer(training_j)
dtm_training$j<-append(dtm_training$j,training_j)
training_v<-as.numeric(training_v)
dtm_training$v<-append(dtm_training$v,training_v)

test_i<-as.integer(test_i)
dtm_test$i<-append(dtm_test$i,test_i)
test_j<-as.integer(test_j)
dtm_test$j<-append(dtm_test$j,test_j)
test_v<-as.numeric(test_v)
dtm_test$v<-append(dtm_test$v,test_v)