require("tm")
require("topicmodels")
require("slam")
i <- read.table("i.txt", header=T, quote="\"")
j <- read.table("j.txt", header=T, quote="\"")
v <- read.table("v.txt", header=T, quote="\"")
time <- read.table("time.txt", header=T, quote="\"")
words <- read.table("words.txt", header=T, quote="\"")
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