

rdirichlet <- function(a) {
    y <- rgamma(length(a), a, 1)
    return(y / sum(y))
}


##################################################################################
#
# Simulate data from topic model, original version
#
#    simulated response is lognormal 
#
##################################################################################
poisson.topic.model <- function() { 
	
	K      <-  100				#	number attributes

	n.doc  <- 2500				#	number documents/listings  
	
	                       
								#	Z is binary matrix that identifies the attributes in each doc
								#	prob controls number of topics in a document (affects mean length)
	Z       <- matrix( rbinom(n.doc * K, size=1, prob=0.3), nrow=n.doc, ncol=K)

	beta    <- rnorm(K)			#	intercept is zero
	sigma   <- 1
	y       <- Z %*% beta + rnorm(n.doc,sd=sigma) 
	summary(regr <- lm(y ~ Z))	#	check "true" model; R2 around 90%


	
	lambda  <- 1+exp(beta)		# say more about positive attributes; scale to get mean length right
	
	
	n.vocab <- 500				# size of vocabulary
								
								# matrix with distributions over vocab for each attribute
	P       <- matrix(0, nrow=K, ncol=n.vocab)
	n.common <-     100 					# number words shared over attributes
	common <- c(rdirichlet( rep(2,n.common) ), rep(0,n.vocab-n.common) ) 
	alpha <- rep(.05, n.vocab)				# dirichlet parm, small alpha makes spiky
	for(k in 1:K) P[k,] <- (common + rdirichlet(alpha))/2
	if(n.common==0) P <- 2 * P
	round(P[1:3,1:15],3); 					# take a peek at three attribute vocabularies
	# apply(P,1,sum)						# prob dist so sum to 1


	W <- matrix(0, nrow=n.doc, ncol=n.vocab)	# word frequencies
	one <- rep(1,n.vocab)
	for(i in 1:n.doc) {
		indx <- which(Z[i,] == 1)
		W[i,] <- rpois(one,lambda[indx] %*% P[indx,])
	}
	
	# check overall lengths
	doc.len <- apply(W,1,sum)
	mean(doc.len); hist(doc.len, breaks=30)
	
	# plots of y on length
	plot(y ~ doc.len); cor(y, doc.len)  		# correlated
	fit <-loess(y ~ doc.len, span=0.3)    # nothing nonlinear
	lines(predict(fit)[order(doc.len)],col="red")   # nothing nonlinear

	lines(lowess(doc.len,y, f=0.3),col="blue")
	
	# check vocab for Zipf distribution ... needs much more skew, more very rare types
	freq <- apply(W,2,sum)
	hist(freq)
	
	o <- order(freq,decreasing=TRUE)	
	W.ordered <- W[,o]
	plot(log(freq[o]) ~ log(1:n.vocab))			# want slope -1
	
	# word regressions (lsa style, with word types in order of overall frequency
	sr <- summary(regr <- lm(y ~ W.ordered[,1:100])); sr
	
	
}




##################################################################################
#
# Simulate data from topic model, original version
#
##################################################################################
original.topic.model <- function() { 

# --- functions

rdirichlet <- function(a) {
    y <- rgamma(length(a), a, 1)
    return(y / sum(y))
}

word.indices <- function(topics, P) { 
	# generate word by sampling distribution of topics
	M <- ncol(P); m<-length(topics) 
	z<-rep(0,m); 
	for(t in 1:m) z[t]<-sample(M,1,prob=P[topics[t],])
	z
}

# --- constants

beta <- c(1, 1, 1, 2, 2, 2, 3, 3,-2,-3) 	# topic weights
K    <- length(beta)		# number of topics

M <- 2000					# word types in vocabulary
n <- 6000					# num of documents

# --- topic distributions over words in K x M matrix P
P <- matrix(0, nrow=K, ncol=M)
                                              # generate words shared by topics
n.common <-     0                             # number shared words 
common <- c(rdirichlet( rep(2,n.common) ), rep(0,M-n.common) ) 
alpha <- rep(.05, M)			                  # dirichlet parms, small alpha are spiky
for(k in 1:K) P[k,] <- (common + rdirichlet(alpha))/2
if(n.common==0) P <- 2 * P
round(P[1:10,1:12],3); apply(P,1,sum)         # prob dist so sum to 1

#     check dependence/correlation among distributions
udv <- svd(P); plot(udv$d); 
#     common words along diagonal
round(cor(t(P)),2);  plot(jitter(P[1,]),jitter(P[2,]), main="Scatter of Two Dist")
round(udv$u[,1],4); 
plot(udv$v[,1])    # common words are the leading n.col values


# --- Z[i,] is distribution of topics in document i (random Dirichlet)
Z <- matrix(0, nrow=n, ncol=K);
for(i in 1:n) Z[i,] <- rdirichlet( rep(1/10,K) )

#     Y is the response
Y <- Z %*% beta + rnorm(n,sd=0.5)

#     check "true" model; sum of the X's = 1 so no intercept, R2 around 90%
summary(regr <- lm(Y ~ Z - 1))


# --- simulate documents; lengths of the documents (neg bin)
vocab <- paste("w",1:M,sep="")
m <- rpois(n,lambda=rgamma(n,30,1))

# words in doc assigned to topic, then pick word for that topic
docs <- rep("",n)
n.topics <- rep(0,n)
for(i in 1:n) {
	topics <- sample(K, m[i], replace=TRUE, prob=Z[i,]); 
	n.topics[i] <- length(unique(topics))
	docs[i] <- paste(round(Y[i],3), paste(vocab[word.indices(topics,P)], collapse=" "), sep=" ")
	}
# count avg number unique topics
hist(n.topics); mean(n.topics)
	
# write words to file with response
write(docs,"/Users/bob/C/text/text_src/sim/sim.txt")
}

