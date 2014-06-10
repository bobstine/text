##################################################################################
#
#  This version uses the Blei-McAuliffe sLDA model
#
# Given α and K topics with distributions β and regression coef η, var σ2
# Generate a document as follows...
# 1. Draw topic proportions θ | α ∼ Dir(α). 
# 2. For each word w_n
# 		(a) Draw topic assignment z_n |θ ∼ Mult(θ)		z_n indicates topic
# 		(b) Draw word w_n|z_n,β_{1:K} ∼ Mult(β_{z_n})		
# 3. Draw response variable 
# 		y|z_{1:N},η,σ2 ∼ N(η'z-bar,σ2).
#  z-bar is the average of the z's for document, the average topic assignment 
#  the topic composition of all words in the document. (proportions of topics)
#  Regressors are empirical rather than expected topic proportions.
#
##################################################################################

source("/Users/bob/C/text/functions.R")

##################################################################################
#
#   Supervised topic model generative process (Blei and McAullife)
#
##################################################################################


par(mfrow=c(1,2))				# [P.pdf] 
	a <- 0.01;
	p1 <- rdirichlet(rep(a,n.vocab)); 	p2 <- rdirichlet(rep(a,n.vocab))
	plot(p1,p2, xlab=expression("P"[1]),ylab=expression("P"[2]), 
				main=expression(paste(alpha,"=0.01")))
	a <- 0.10;
	p1 <- rdirichlet(rep(a,n.vocab)); 	p2 <- rdirichlet(rep(a,n.vocab))
	plot(p1,p2, xlab=expression("P"[1]),ylab=expression("P"[2]), 
				main=expression(paste(alpha,"=0.10")))
reset()
	
	n.vocab <- 1500
	K <- 30								# number of topics
	P <-matrix(0,nrow=K,ncol=n.vocab)	# dist over words for each topic
	
	alpha.P <- 0.10
	for(i in 1:K) P[i,] <- rdirichlet(rep(alpha.P,n.vocab))	
	
	plot(sqrt(P[1,]),sqrt(P[2,]))     	# disjoint if alpha = 0.01, more common if .1
	plot(P[3,])
	
	alpha <- rep(0.4,K)					# mix of topics within documents
	n <- 5000							# documents
	avg.len <- 50						# avg document length
	theta <- matrix(0,nrow=n,ncol=K)	# expected topic mix
	Z	  <- matrix(0,nrow=n,ncol=K)	# observed topic mix
	doc.len <- rep(0,n)
	for(i in 1:n) {						
		doc.len[i] <- rpois(1,avg.len)
		theta[i,] <- rdirichlet(alpha)		
		Z[i,] <- as.vector(rmultinom(1,doc.len[i],theta[i,]))
	}
	
	W <- matrix(0,nrow=n, ncol=n.vocab)	# generate y and W
	eta <- rnorm(K, mean=2, sd=1) 		# topic regr coefs, need some neg
	mu.y <- rep(0,n)
	for(i in 1:n) {						
		for(k in 1:K) W[i,]<-W[i,]+rmultinom(1,Z[i,k],P[k,])
		mu.y[i] <- sum(eta * Z[i,])
	}
	
	target.r2 <- 0.6
	Y <- mu.y + rnorm(n, sd=sqrt((1-target.r2) * var(mu.y) / target.r2))
	summary(regr.z <-	lm(Y ~ Z) )
	
	
# --- Zipf?

	freq <- apply(W,2,sum)
	sort.freq <- sort(freq[freq>0],decreasing=TRUE) 
	lf <- log(sort.freq); lr <- log(1:length(sort.freq))
	plot(sort.freq,log="xy", xlab="Word Rank",ylab="Frequency")
	lf <- lf[1:100]; lr <- lr[1:100]
	regr <- lm(lf ~ lr); coefficients(summary(regr)); 
	lines(exp(predict(regr)), col="red")
	cat(sum(freq==0)," unused words\n");

	
##################################################################################
#
#   Fit regression models
#
##################################################################################

# --- fit using expected(theta) and actual(z) topic proportions
	hist(Y)
	
	# --- these are independent rather than correlated if use z-bar
	#	  but become very correlated if generate with Z if all eta>0
	#	  so use some negative etas, but then cannot regress on topic mix
	#	  and have to weight up by doc length
	
	plot(Y ~ doc.len); cor(Y,doc.len)^2
	
	summary(rt <- lm(Y ~ theta) )
	summary(rz <- lm(Y ~ Z) );

# --- ordered words, drop those absent
	freq <- apply(W,2,sum)
	o <- order(freq,decreasing=TRUE); cat("Most common words: ", o[1:20], "\n");
	o <- o[freq[o]>0]
	W.ordered <- W[,o]
	dim(W.ordered)
	
# --- word regressions
	# saturated model to get max possible R2 with words + len
	sr <- summary(regr.W <- lm(Y ~ W.ordered[,1:200]));  sr
	coef.summary.plot(sr, "Word Frequencies", omit=2)

	predictive.r2(rz)
	predictive.r2(regr.W)

# --- SVD 				  **** Be sure to generate W.ordered first  ****

	# --- LSA analysis, raw counts
	udv <- svd(W.ordered)
	plot(udv$d[1:100], log="xy")
	U <- udv$u	
	
	# --- LSA analysis, CCA scaled
	W.cca <- (1/sqrt(doc.len)) * W.ordered
	w.freq <- apply(W.ordered,2,sum)
	W.cca <- t( (1/sqrt(w.freq)) * t(W.cca))
	udv.cca <- svd(W.cca)
	plot(udv.cca$d[1:100], log="xy")
	U.cca <- udv.cca$u	
	
par(mfrow=c(1,2))				# [spectra.pdf]
	plot(udv$d[1:100], log="xy", xlab="Component", ylab="Singular Value",
		main="Raw Frequencies")
	plot(udv.cca$d[1:100], log="xy", xlab="Component", ylab="Singular Value",
		main="CCA Normalization")
reset()

# --- leading singular value determined by word frequency
f  <- freq[o]; f <- f[f>0]
ni <- rowSums(W);

plot(f, udv$v[,1])
plot(f, udv.cca$v[,1])
plot(ni, udv$u[,1])
		
# --- SVD regressions
	sr <- summary(regr.u <- lm(Y ~ U[,1:100]));  sr
	coef.summary.plot(sr, "Singular Vectors U", omit=2)
	
	sr <- summary(regr.cca <- lm(Y ~ U.cca[,1:100]));  sr
	coef.summary.plot(sr, "Singular Vectors U, cca", omit=2)

	r2<- rbind(predictive.r2(rz), predictive.r2(regr.cca), 
	               predictive.r2(regr.u), predictive.r2(regr.W))
	r2
	round(r2[,3],3)



##################################################################################
#
# Simulate topic model for paper  (package above into functions)
#
##################################################################################

# --- globals

n.vocab <- 1500	
avg.len <- 50				# avg document length

K <- 30						# number of topics
n <- 5000					# number of documents

# --- functions

gen.data <- function(alpha.P) {
	P <-matrix(0,nrow=K,ncol=n.vocab)	# dist over words for each topic
	for(i in 1:K) P[i,] <- rdirichlet(rep(alpha.P,n.vocab))	
	alpha <- rep(0.4,K)					# mix of topics within documents
	Z	  <- matrix(0,nrow=n,ncol=K)	# observed topic mix
	doc.len <- rpois(n,avg.len)
	for(i in 1:n) {	Z[i,] <- as.vector(rmultinom(1,doc.len[i],rdirichlet(alpha))) }
	W <- matrix(0,nrow=n, ncol=n.vocab)	# generate y and W
	eta <- rnorm(K, mean=2, sd=1) 		# topic regr coefs, need some neg
	mu.y <- rep(0,n)                   
	for(i in 1:n) {						
		for(k in 1:K) W[i,]<-W[i,]+rmultinom(1,Z[i,k],P[k,])
		mu.y[i] <- sum(eta * Z[i,])
	}
	target.r2 <- 0.6
	Y <- mu.y + rnorm(n, sd=sqrt(var(mu.y) * (1-target.r2)/target.r2))
	freq <- apply(W,2,sum)
	o <- order(freq,decreasing=TRUE)
	o <- o[freq[o]>0]
	W.ordered <- W[,o]
	list(W=W.ordered,Y=Y, Z=Z)
}

fit.regressions <- function(data, k.W, k.U) {
	regr.Z <- lm(data$Y ~ data$Z)
	regr.W <- lm(data$Y ~ data$W[,1:k.W])
	udv <- svd(data$W)
	regr.u <- lm(data$Y ~ udv$u[,1:k.U])
	doc.len <- rowSums(data$W)
	W.cca <- (1/sqrt(doc.len)) * data$W
	w.freq <- colSums(data$W)
	W.cca <- t( (1/sqrt(w.freq)) * t(W.cca))
	udv <- svd(W.cca)
	regr.cca <- lm(data$Y ~ udv$u[,1:k.U])
	(rbind(predictive.r2(regr.Z), predictive.r2(regr.cca), 
	      predictive.r2(regr.u), predictive.r2(regr.W)))[,3]
}

doit <- function(filename, n.reps, alpha) {
	# write predictive R2s to file, one line at a time
	for(rep in 1:n.reps) {
		data <- gen.data(alpha);
		results <- fit.regressions(data, 200, 100)
		cat(rep," ")
		write(results, filename, append=TRUE);
	}
}

system.time( doit("/Users/bob/Desktop/sLDA_01.txt", 100, 0.01) )

system.time( doit("/Users/bob/Desktop/sLDA_10.txt", 100, 0.10) )

#### ---- Plot simulation results
R2.10 <- read.table("~/Desktop/sLDA_10.txt")

plot(R2.10[,1], type="l", ylim=c(0,0.6))
lines(R2.10[,2])
lines(R2.10[,3])
lines(R2.10[,4])

R2.01 <- read.table("~/Desktop/sLDA_01.txt")

plot(R2.01[,1], type="l", ylim=c(0.35,0.65))
lines(R2.01[,2])
lines(R2.01[,3])
lines(R2.01[,4])

pairs(R2.01)

