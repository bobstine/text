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

source("~/C/text/functions.R")


##################################################################################
#
#   Supervised topic model generative process (Blei and McAullife)
#
#		P holds topic distributions over words (rows)
#		Sort P and Z (observed counts) to have largest first
#
##################################################################################
	
	n.vocab <- 1500
	K <- 30								# number of topics
	P <-matrix(0,nrow=K,ncol=n.vocab)	# dist over words for each topic
	
	alpha.P <- 0.01   					# smaller alpha implies more diffuse, vary over topics
	for(i in 1:K) P[i,] <- rdirichlet(rep(alpha.P,n.vocab))	# small alpha implies highly skewed
	P <- P[,order(colSums(P), decreasing=TRUE)]				# sort so big prob are first
	
	plot(sqrt(P[1,]),sqrt(P[2,]))     	# disjoint if alpha = 0.01, more common if .1
	plot(P[3,])							# weights on specific words

										# plot of example topic distributions
	par(mfrow=c(1,2))					# [P.pdf] 
		a <- 0.01;
		p1 <- rdirichlet(rep(a,n.vocab)); 	p2 <- rdirichlet(rep(a,n.vocab))
		plot(p1,p2, xlab=expression("P"[1]),ylab=expression("P"[2]), 
				main=expression(paste(alpha,"=0.01")))
		a <- 0.10;
		p1 <- rdirichlet(rep(a,n.vocab)); 	p2 <- rdirichlet(rep(a,n.vocab))
		plot(p1,p2, xlab=expression("P"[1]),ylab=expression("P"[2]), 
				main=expression(paste(alpha,"=0.10")))
	reset()

	alpha <- rep(0.4,K)					# mix of topics within documents
	n <- 5000							# documents
	avg.len <- 50						# avg document length
	theta <- matrix(0,nrow=n,ncol=K)	# expected topic mix
	Z	  <- matrix(0,nrow=n,ncol=K)	# number of words from each topic
	doc.len <- rep(0,n)
	for(i in 1:n) {						
		doc.len[i] <- rpois(1,avg.len)
		theta[i,] <- rdirichlet(alpha)		
		Z[i,] <- as.vector(rmultinom(1,doc.len[i],theta[i,]))
	}
	Z <- Z[order(rowSums(Z), decreasing=TRUE),]
	
	W <- matrix(0,nrow=n, ncol=n.vocab)	# generate y and W, and 
	W.ev <- W;							#  expected value of W
	eta <- rnorm(K, mean=2, sd=1) 		# topic regr coefs, need some neg
	mu.y <- rep(0,n)
	for(i in 1:n) {	
		W.ev[i,] <- doc.len[i] * theta[i,] %*% P # prob distribution over words for each doc		
		for(k in 1:K) W[i,]<-W[i,]+rmultinom(1,Z[i,k],P[k,])
		mu.y[i] <- sum(eta * Z[i,])		# words within a topic are exchangeable
	}
	
	
	# --- check that W counts are centered on their means for doc with lots
	plot(rowSums(W.ev), rowSums(W))		# perfect corr since use Poisson numbers
	
	i <- 1; plot(W.ev[i,], W[i,], xlab="Expected Word Counts", ylab="Observed")
	points( W.ev[i,], fitted.values(lm(W[i,] ~ poly(W.ev[i,],8)) ), col="gray")
	
	
# --- Zipf?  Right initial shape, but not nearly so steep as needed (power ≈ -0.3)

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
#   Analyze term/document counts from topic model
#
##################################################################################

	# --- LSA analysis, raw counts
	udv <- svd(W)
	U <- udv$u
	V <- udv$v
	
	# --- LSA analysis, expected counts
	#		singular values from ev much more distinctive
	udv.ev <- svd(W.ev)
	U.ev <- udv.ev$u
	V.ev <- udv.ev$v
	
	par(mfrow=c(2,1))
		plot(udv   $d[1:100], log="xy", ylab="Singular Values")
		plot(udv.ev$d[1:100], log="xy", ylab="Singular Values")
	reset()
	
	# --- LSA analysis, CCA scaled
	W.cca   <- (1/sqrt(ni)) * W
	W.cca   <- t( (1/sqrt(fj)) * t(W.cca) )
	udv.cca <- svd(W.cca)
	plot(udv.cca$d[1:100], log="xy")
	U.cca <- udv.cca$u	
	V.cca <- udv.cca$v
	
	# --- Square root of counts (stabilize Poisson)
	udv.sr <- svd(sqrt(W))
	plot(udv.sr$d[1:100], log="xy")
	U.sr <- udv.sr$u
	V.sr <- udv.sr$v
	
	par(mfrow=c(2,1))				# [spectra.pdf]
		plot(udv$d[1:100], log="xy", xlab="Component", ylab="Singular Value",
			main="Raw Frequencies")
		plot(udv.cca$d[1:100], log="xy", xlab="Component", ylab="Singular Value",
			main="CCA Normalization")
	reset()

# --- leading singular vectors determined by word frequency
	par(mfrow=c(1,2))
		plot(fj, udv.ev $v[,1])
		plot(ni, udv.ev $u[,1])
	reset()
	par(mfrow=c(1,2))
		plot(fj, udv    $v[,1])
		plot(ni, udv    $u[,1])
	reset()
	par(mfrow=c(1,2))
		plot(fj, udv.cca$v[,1])
		plot(ni, udv.cca$u[,1])
	reset()
	par(mfrow=c(1,2))
		plot(fj, udv.sr $v[,1])
		plot(ni, udv.sr $u[,1])
	reset()

# --- structure of singular vectors, loadings highlight those with large magnitude overall
#		cca components have less 'pointy' loadings
	j <- 4; k <- j+1;  
	
	plot(V.ev[,j],V.ev[,k], xlab=paste0("V_ev(",j,")"), ylab=paste("V_ev(",k,")"))
	
	plot(V[,j],V[,k], xlab=paste0("V_raw(",j,")"), ylab=paste("V_raw(",k,")"))
	
	plot(V.sr[,j],V.sr[,k], xlab=paste0("V_sr(",j,")"), ylab=paste("V_sr(",k,")"))

	plot(V.cca[,j],V.cca[,k], xlab=paste0("V_cca(",j,")"), ylab=paste("V_cca(",k,")"))


# --- relationship between loading spaces
#		expected values
	cc   <- cancor(V    [,1:50], V.ev[,1:50]); plot(cc$cor)
	cc.1 <- cancor(V.cca[,1:50], V.ev[,1:50]); points(cc.1$cor, col="red", cex=.5)
	cc.2 <- cancor(V.sr [,1:50], V.ev[,1:50]); points(cc.2$cor, col="blue", cex=.5)
#		topic distributions
	cc   <- cancor(V    [,1:50], P[,1:50]); plot(cc$cor)
	cc.1 <- cancor(V.cca[,1:50], V.ev[,1:50]); points(cc.1$cor, col="red", cex=.5)
	cc.2 <- cancor(V.sr [,1:50], V.ev[,1:50]); points(cc.2$cor, col="blue", cex=.5)

##################################################################################
#
#   Regression models
#
##################################################################################

	
# --- check that W regression is in the right ballpark
	target.r2 <- 0.6
	Y <- mu.y + rnorm(n, sd=sqrt((1-target.r2) * var(mu.y) / target.r2))
	summary(regr.z <-	lm(Y ~ Z) )
	
# --- word regressions
	# saturated model to get max possible R2 with words + len
	sr <- summary(regr.W <- lm(Y ~ W[,1:200]));  sr
	coef.summary.plot(sr, "Word Frequencies", omit=2)

	predictive.r2(rz)
	predictive.r2(regr.W)

# --- fit using expected(theta) and actual(z) topic proportions
	hist(Y)
	
	# --- these are independent rather than correlated if use z-bar
	#	  but become very correlated if generate with Z if all eta>0
	#	  so use some negative etas, but then cannot regress on topic mix
	#	  and have to weight up by doc length
	
	plot(Y ~ doc.len); cor(Y,doc.len)^2
	
	summary(rt <- lm(Y ~ theta) )
	summary(rz <- lm(Y ~ Z) );

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


