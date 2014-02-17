
source("/Users/bob/C/text/functions.R")

rdirichlet <- function(a) {
    y <- rgamma(length(a), a, 1)
    return(y / sum(y))
}


##################################################################################
#
# Simulate data from topic model, Poisson version
#
#    simulated response is lognormal (eg, log price)
#
##################################################################################
poisson.topic.model <- function() { 
	
	n.doc  <- 2500					#	number observed documents/listings  
	K      <-   50					#	number possible attributes of varying expected values
									
	x <- sort(rexp(K))				#	value of attributes, lognormal, increasing mu
	mu <- x^2			
	hist(mu, main="Simulated Contributions to Log Prices", breaks=50)		
	
									#   number attributes in each document is 1 + neg binomial
	A.sum   <- 1 + rpois(rep(1,n.doc), 4*rgamma(n.doc,shape=2))
	hist(A.sum); mean(A.sum)		#	average number of attributes
									#	A identifies attributes in each doc	
	A <- matrix(0, nrow=n.doc, ncol=K)
	for(i in 1:n.doc) A[i,] <- rmultinom(1, A.sum[i], rep(1,K))

	Y <- A %*% mu					#	perfect true model
	hist(Y, main="Dist of Log Price")
	plot(Y ~ A.sum, xlab="Number Latent Attributes", ylab="Log Price"); cor(Y,A.sum)^2

	n.vocab <- 1500					# matrix with distributions over vocab for each attribute
	P       <- matrix(0, nrow=K, ncol=n.vocab)
	for(k in 1:K) P[k,] <- rdirichlet(rep(.1, n.vocab))
	round(P[1:3,1:15],3); 			# take a peek at three attribute vocabularies
	apply(P,1,sum)[1:4]				# prob dist so sum to 1
	plot(P[1,],P[2,])

	# --- generate doc/word matrix
	W <- matrix(0, nrow=n.doc, ncol=n.vocab)	#	word frequencies
	lambda <- 6									#	overall expected words per attribute; random effect
	one <- rep(1,n.vocab)
	m <- A %*% P;
	for(i in 1:n.doc) { W[i,] <- rpois(one,rgamma(n.vocab, lambda/5, scale=5)*m[i,]) }
	cat("Max word frequencies in first 5 docs", apply(W,1,max)[1:5])
	doc.len <- apply(W,1,sum)					# check document lengths
	mean(doc.len); hist(doc.len, breaks=30)
	
	# --- check vocab for Zipf distribution
	freq <- apply(W,2,sum)
	sort.freq <- sort(freq[freq>0],decreasing=TRUE) 
	lf <- log(sort.freq); lr <- log(1:length(sort.freq))
	plot(lf ~ lr, xlab="log rank",ylab="log frequency")		# want slope -1
	lf <- lf[1:300]; lr <- lr[1:300]
	regr <- lm(lf ~ lr); summary(regr); lines(lr, predict(regr), col="red")
	cat(sum(freq==0)," unused words\n");
	
	# --- regress log price on length
	plot(Y ~ doc.len, xlab="Doc Length", ylab="Log Price");  
	fit <-loess(Y ~ doc.len, span=0.3)    			# nothing nonlinear
	o <- order(doc.len)
	lines(doc.len[o], predict(fit)[o],col="red")	# nothing nonlinear
	r <- lm(Y ~ doc.len); cor(Y, doc.len)^2
	lines(doc.len[o], predict(r)[o], col="blue")
	
		
	# word regressions (word types in order of overall frequency)
	W.ordered <- W[,order(freq,decreasing=TRUE)]
	sr <- summary(regr <- lm(Y ~ W.ordered[,1:100]));  sr
	
	# plot regression coefs
	par(mfrow=c(1,2))  
		y <- abs(sr$coefficients[-1,3])
		x <- 1:length(y)        # some may be singular
		threshold <- -qnorm(.025/length(y))
		plot(x,y,	xlab= "Columns of W", ylab="|t|", main="")
		abline(h=threshold, col="gray", lty=4)
		abline(h=sqrt(2/pi), col="cyan")
		lines(lowess(x,y,f=0.3), col="red")
		half.normal.plot(y)
	reset()
	
# LSA analysis

	udv <- svd(W[,freq>0])
	U <- udv$u
	
	plot(udv$d[1:300], log="y")
	
	lsa.sr <- summary(regr <- lm(Y ~ U[,1:100])); lsa.sr
	
	# plot regression coefs
	par(mfrow=c(1,2))  
		y <- abs(lsa.sr$coefficients[-(1:2),3])
		x <- 1:length(y)        # some may be singular
		threshold <- -qnorm(.025/length(y))
		plot(x,y,	xlab= "LSA Variables", ylab="|t|", main="")
		abline(h=threshold, col="gray", lty=4)
		abline(h=sqrt(2/pi), col="cyan")
		lines(lowess(x,y,f=0.3), col="red")
		half.normal.plot(y)
	reset()
	

	
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

