#
#  Analysis of ANES data from LDC
#

dataPath <- "/Users/bob/data/text/text_src/anes/"
add.path <- function(file) { paste(dataPath,file,sep="") }

#####################################################################################
#
#    Office recognition for Gordon Brown
#
#####################################################################################

Brown <- read.csv(add.path("brown.csv")); dim(Brown)   # 2098 x 10, but last was empty

# --- text is in second column; convert to lower-case text

typeof(Brown$InterviewerTranscript)

Brown$InterviewerTranscript <- tolower(as.character(Brown$InterviewerTranscript))
typeof(Brown$InterviewerTranscript)


# --- separate into words

Brown$words <- strsplit(Brown$InterviewerTranscript," ")
Brown$words[[1]]

Brown$nTokens <- as.vector( sapply(Brown$words,length) )
Brown$words[[which.max(Brown$nTokens)]]

hist(Brown$nTokens)
cat("Median number of tokens", median(Brown$nTokens))


# --- build vocabulary

all.counts <- table(c(Brown$words, recursive=TRUE))
all.counts <- sort(all.counts, decreasing=TRUE)
all.counts[1:20]
all.types <- names(all.counts)


# --- remove rare types (oov); 680 for brown

oov.types <- all.types[(all.counts==1)|(all.counts == 2)]
length(oov.types)

types <- all.types[all.counts>2]; length(types)


# --- create document/type matrix W

equal <- function(a,b) { as.numeric(a == b) }

#     test it
i <- 2;
Brown$InterviewerTranscript[[i]]
w <- (outer(Brown$words[[i]], types, FUN="equal"))
b <- apply(w,2,sum)
sum(b)

W <- matrix(0, nrow=nrow(Brown), ncol=length(types))
for(i in 1:nrow(Brown)) W[i,] <- apply(outer(as.vector(Brown$words[[i]]),types,FUN="=="),2,sum)

#    difference in these counts are oov tokens
apply(W,1,sum)[1:10]
Brown$nTokens[1:10]

#    tack on oov
W <- cbind(W,Brown$nTokens-apply(W,1,sum))
colnames(W)<-c(types,"OOV")

#    check that lengths match total counts 10453
sum(W)
sum(Brown$nTokens)


# --- LSA , after frequency and entropy normalization

entropy <- function(p) {    # 0 to 1
	p.norm <- p[p>0]/sum(p)
	-sum(log(p.norm)*p.norm)/log(length(p)) }
	
	
L <- W;
for(i in 1:nrow(L)) L[i,] <- L[i,]/sum(L[i,])
for(j in 1:ncol(L)) L[,j] <- (1-entropy(L[,j])) * L[,j]

# find missing for(i in 1:nrow(L)) if(any(is.na(L[i,]))) cat(i,"\n")

udv <- svd(L)

plot(udv$d[-length(udv$d)], log="y",
	main="Singular Values, LSA (Brown Office Recognition)", xlab="Index", ylab="Singular Value")


# --- Cluster documents using kmeans (0,10=other; 1,11=dk)

k <- 2
clusters <- kmeans(udv$u[,1:100], centers = k)
clusters$size
100 * clusters$betweenss/clusters$totss

table(clusters$cluster,Brown$Code0)

Brown$InterviewerTranscript[clusters$cluster==2]
Brown[clusters$cluster==2,c(2,3)]

