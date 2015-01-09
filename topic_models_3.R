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
	
	n.outlier.words <- 20				# reserved for words only in outlier listing (optional)
	n.vocab <- 1000 					# total vocab
	K <- 30								# number of topics
	P <-matrix(0,nrow=K,ncol=n.vocab-n.outlier.words)	# dist over words for each topic; pad with 0
	
# --- generate topic distributions
	alpha.P <- 0.02   					# smaller alpha implies high skew, less overlap
	for(i in 1:K) P[i,] <- rdirichlet(rep(alpha.P,n.vocab-n.outlier.words))	
	P <- P[,order(colSums(P), decreasing=TRUE)]				# sort so big total prob are first
	P <- cbind(P, matrix(0,nrow(P), n.outlier.words))

	plot(P[3,])							# weights on specific words
	plot(sqrt(P[1,]),sqrt(P[2,]))     	# disjoint if alpha = 0.01, some common if .1

	par(mfrow=c(1,2))					# plot of example topic distributions [P.pdf] 
		a <- 0.05; 
		p1 <- rdirichlet(rep(a,n.vocab)); 	p2 <- rdirichlet(rep(a,n.vocab))
		plot(p1,p2, xlab=expression("P"[1]),ylab=expression("P"[2]), 
				main=expression(paste(alpha,"=0.05")))
		a <- 0.10;
		p1 <- rdirichlet(rep(a,n.vocab)); 	p2 <- rdirichlet(rep(a,n.vocab))
		plot(p1,p2, xlab=expression("P"[1]),ylab=expression("P"[2]), 
				main=expression(paste(alpha,"=0.10")))
	reset()

# --- generate samples from mixture of topics
	alpha <- rep(0.4,K)					# mix of topics within documents
	n <- 5000							# documents
	avg.len <- 50						# avg document length
	theta <- matrix(0,nrow=n,ncol=K)	# expected topic mix
	Z	  <- matrix(0,nrow=n,ncol=K)	# number of words from each topic
	doc.len <- sort(rpois(n,avg.len), decreasing=TRUE) 
	for(i in 1:n) {						
		theta[i,] <- rdirichlet(alpha)		
		Z[i,] <- as.vector(rmultinom(1,doc.len[i],theta[i,]))
	}
	plot(theta[1,], xlab="Topic", ylab="Share of Vocabulary", main="Topic Mix")
	
	W <- matrix(0,nrow=n, ncol=n.vocab)	# generate y and W, and 
	W.ev <- W;							#  expected value of W
	eta <- rnorm(K, mean=2, sd=1) 		# topic regr coefs, need some neg
	mu.y <- rep(0,n)
	for(i in 1:n) {	
		W.ev[i,] <- doc.len[i] * theta[i,] %*% P # prob distribution over words for each doc		
		for(k in 1:K) if(Z[i,k]>0) W[i,]<-W[i,]+rmultinom(1,Z[i,k],P[k,])
		mu.y[i] <- sum(eta * Z[i,])		# words within a topic are exchangeable
	}
	
	# --- check that W counts are centered on their means for doc with lots
	plot(rowSums(W.ev), rowSums(W))		# perfect corr since use Poisson numbers
	
	# --- within document regr of counts on ev
	#		p-values are not correct since counts are not independent?
	i <- 2; 
	summary(r1 <- glm(W[i,] ~      W.ev[i,]   , family="poisson"))
	summary(rp <- glm(W[i,] ~ poly(W.ev[i,],5), family="poisson"))
	plot(W.ev[i,], W[i,], xlab="Expected Word Counts", ylab="Observed")
	abline(a=0,b=1,col="gray")
	lines(lowess(W.ev[i,], W[i,],f=0.2),col="cyan")
	points( W.ev[i,], fitted.values(r1), col="red")
	points( W.ev[i,], fitted.values(rp), col="gray")	
	
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

	# --- marginals  Leave zeros for outlier insertion/handling
	fj <- colSums(W)
	ni <- rowSums(W)  # = doc.len
	
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
		plot(udv.ev$d[1:100], log="xy", ylab="Singular Values, EV")
	reset()
	
	# --- LSA analysis, CCA scaled
	#		gap is evident... but this will be problem later
	recip.sqrt <- function(n) { if(n>0) return(1/sqrt(n)) else return(1)}
	W.cca   <- sapply(ni,recip.sqrt) * W
	W.cca   <- t( sapply(fj,recip.sqrt) * t(W.cca) )
	udv.cca <- svd(W.cca)
	plot(udv.cca$d[1:100], log="xy")
	U.cca <- udv.cca$u	
	V.cca <- udv.cca$v
	
	# --- Square root of counts (stabilize Poisson)
	#		gap is pretty hard to see
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

# --- leading singular vectors u_1,v_1 determined by word frequency
	par(mfrow=c(1,2))
		plot(fj, udv.ev $v[,1])
		plot(ni, udv.ev $u[,1])
	reset()
	par(mfrow=c(1,2))
		plot(fj, udv    $v[,1], xlab=expression("f"["j"]),ylab=expression("V"[1]))
		plot(ni, udv    $u[,1], xlab=expression("n"["i"]),ylab=expression("U"[1]))
	reset()
	# CCA produces 'mathematical' fit
	par(mfrow=c(1,2))
		plot(sqrt(fj),-udv.cca$v[,1], xlab=expression(sqrt("f"["j"])),ylab=expression("V"[1]))
		plot(sqrt(ni),-udv.cca$u[,1], xlab=expression(sqrt("n"["i"])),ylab=expression("U"[1]))
	reset()
	udv.cca$d[1:5];
	v1 <- udv.cca$v[,1]; summary( lm(v1 ~ sqrt(fj)) )
	u1 <- udv.cca$u[,1]; summary( lm(u1 ~ sqrt(ni)) )
	#
	par(mfrow=c(1,2))
		plot(fj, udv.sr $v[,1])
		plot(ni, udv.sr $u[,1])
	reset()

# --- structure of singular vectors, loadings highlight those with large magnitude overall
#		cca components have less 'pointy' loadings
	j <- 4; k <- j+1;  
	
	plot(V.ev[,j],V.ev[,k],   xlab=paste0("V_ev(",j,")"), ylab=paste("V_ev(",k,")"))
	
	plot(V[,j],V[,k],         xlab=paste0("V_raw(",j,")"), ylab=paste("V_raw(",k,")"))
	
	plot(V.sr[,j],V.sr[,k],   xlab=paste0("V_sr(",j,")"), ylab=paste("V_sr(",k,")"))

	plot(V.cca[,j],V.cca[,k], xlab=paste0("V_cca(",j,")"), ylab=paste("V_cca(",k,")"))


# --- relationship between loading spaces
#		expected values
	cc   <- cancor(V    [,1:50], V.ev[,1:50]); plot(cc$cor)
	cc.1 <- cancor(V.cca[,1:50], V.ev[,1:50]); points(cc.1$cor, col="red", cex=.5)
	cc.2 <- cancor(V.sr [,1:50], V.ev[,1:50]); points(cc.2$cor, col="blue", cex=.5)
#		topic distributions
	cc   <- cancor(V    [,1:50], t(P[1:10,])); plot(cc$cor)
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
#   Outlier effects on components
#
#						prob want a smallish vocab to run faster
#
##################################################################################

# --- outlier distribution optionally has disjoint words; check L2 dist from others
	P.out <- c(1.00*rdirichlet(rep(alpha.P,n.vocab-n.outlier.words)), 	
			   0.00*rdirichlet(rep(2*alpha.P,n.outlier.words)))
	# R2 < 0.01 with alpha.P = 0.02; closer to 5% with alpha.P = 0.05
	summary(lm(P.out ~ t(P)))  

# --- add extra first row for the outlier  (so can put other outliers without changing rest)
	W <- rbind(rep(0,ncol(W)),W)
	
	W[1,] <- as.vector(rmultinom(1,2*avg.len,P.out))  # add cases to make easier to see?
	ni    <- rowSums(W)
	fj    <- colSums(W)
	
# --- LSA analysis, raw counts
	udv.out <- svd(W)
	U.out <- udv.out$u; V.out <- udv.out$v
	plot  (udv.out$d[1:100], log="xy", xlab="Component", ylab="Singular Value",
			main="Raw Singular Values, Added Topic Outlier", col="blue")
	points(udv    $d[1:100], col="gray", cex=0.5)
	
# --- LSA analysis, CCA scaled with outlier shows another topic before gap
	W.cca.out   <- 1/sqrt(ni) * W
	W.cca.out   <- t( sapply(fj,recip.sqrt) * t(W.cca.out) )
	udv.cca.out <- svd(W.cca.out)
	U.cca.out <- udv.cca.out$u; V.cca.out <- udv.cca.out$v
	plot  (udv.cca.out$d[1:100], log="xy", xlab="Component", ylab="Singular Value",
			main="CCA Normalization, Added Topic", col="blue")
	points(udv.cca    $d[1:100], col="gray", cex=0.5)
	
	# --- svd fits this outlier row separately (regression)
	plot(U.cca.out[,2])  						# all weight on this one document
	cbind(udv.cca.out$d[1:10],udv.cca$d[1:10])	# shifts down by 1
	pairs(cbind(W[1,], -V.cca.out[,2], P.out), c("Data","V","P.out"))
	#		can interpret as a fit?  Yep... So if you can fit this without messing
	#		up the rest, esp if the row offers a particularly large SS
	#		CCA scaling makes for a large SS since fj is low for its words.
	y   <- W.cca.out[1,]; c(sqrt(sum(y^2)), udv.cca.out$d[2])  # sing value is sqrt(ss)
	boxplot(sqrt(rowSums(W.cca.out^2)),log="y")
	fit <- udv.cca.out$d[2] * U.cca.out[1,2] * V.cca.out[,2]
	plot(fit,y); abline(a=0,b=1,col="gray")
	#		ss after scaling
	d <- (W[1,]/sqrt(n[1])) * sapply(fj,recip.sqrt); d <- d - mean(d); sqrt(d %*% d)
	
	# --- How close are LSA and CCA subspaces (before/after outlier)?
	plot( 	cancor(U[,1:75], U.cca[,1:50])$cor   )
	p <- rep(1,50); p[31]<-19
	points( cancor(U.out[,1:75], U.cca.out[,1:50])$cor ,pch=p, col="red"  )
	
	# --- find the topic? highlight words with largest 20 probs in P.out, then CCA for recovery
	color <- rep("gray", ncol(W))
	color[ (order(P.out,decreasing=TRUE))[1:20] ] <- "blue"	
	plot(P.out, V.cca.out[,2], col=color)
	
	cc <- cancor(V.cca.out[,1:50], cbind(t(P),P.out))
	# --- the component that is least well explained is the new one
	plot(cc$cor)
	j<-31; 
	par(mfrow=c(2,1)); 
		plot(cc$xcoef[,j]); abline(h=0,col="gray"); 
		plot(cc$ycoef[,j]); abline(h=0,col="gray");
	reset()
	
	# --- structure of singular vectors
	j <- 2;  
	plot(V.cca.out[,j],V.cca[,j], xlab=paste0("V_out(",j,")"), ylab=paste("V_cca(",j,")"), col=color)
	
	par(mfrow=c(1,2))
		j <- 2; k <- j+1
		plot(V.cca[,j],V.cca[,k], xlab=paste0("V_cca(",j,")"), ylab=paste("V_cca(",k,")"),
			col=color)
		sj <- sign(cor(V.cca[,j],V.cca.out[,j]))
		sk <- sign(cor(V.cca[,k],V.cca.out[,k]))
		plot(sj*V.cca.out[,j],sk*V.cca.out[,k], xlab=paste0("V_out(",j,")"), 
												ylab=paste("V_out(",k,")"),col=color)
	reset()
	

	
# --- now slip in the new words as well and look at spectrum
	P.out.new <- c(0.90*P.out[1:((ncol(W)-n.outlier.words))],	
			       0.10*rdirichlet(rep(2*alpha.P,n.outlier.words)))
	plot(P.out,P.out.new)
	W[n,] <- w.out.new <- as.vector(rmultinom(1,ni[n],P.out.new))
	fj    <- colSums(W)
	plot(W[n,])

	# --- LSA analysis, CCA scaled with outlier shows another topic before gap
	W.cca.out.new <- 1/sqrt(ni) * W
	W.cca.out.new <- t( sapply(fj,recip.sqrt) * t(W.cca.out.new) )
	udv.cca.out.new <- svd(W.cca.out.new)
	plot  (udv.cca.out.new$d[1:100], log="xy", xlab="Component", ylab="Singular Value",
			main="CCA Normalization, Added Topic with New Words", col="red")
	points(udv.cca.out    $d[1:100], col="blue", cex=0.5)
	points(udv.cca        $d[1:100], col="black", cex=0.5)
	U.cca.out.new <- udv.cca.out.new$u	
	V.cca.out.new <- udv.cca.out.new$v


# --- find the topic? highlight words with largest 20 P.out
	color <- rep("gray", ncol(W))
	color[ (order(P.out.new,decreasing=TRUE))[1:20] ] <- "red"	
	plot(P.out.new, V.cca.out.new[,2], col=color); abline(a=0,b=1,col="gray")

	# --- structure of singular vectors, highlight words with large P.out
	color <- rep("black", ncol(W))
	color[ (order(P.out.new,decreasing=TRUE))[1:20] ] <- "red"
	j <- 2;  
	plot(V.cca.out.new[,j],V.cca[,j], xlab=paste0("V_out+new(",j,")"), 
											ylab=paste("V_cca(",j,")"), col=color)
	
	par(mfrow=c(1,2))
		j <- 2; k <- j+1
		plot(V.cca[,j],V.cca[,k], xlab=paste0("V_cca(",j,")"), ylab=paste("V_cca(",k,")"),
			col=color)
		sj <- sign(cor(V.cca[,j],V.cca.out.new[,j]))
		sk <- sign(cor(V.cca[,k],V.cca.out.new[,k]))
		plot(sj*V.cca.out[,j],sk*V.cca.out[,k], xlab=paste0("V_out+new(",j,")"), 
										ylab=paste("V_out+new(",k,")"), col=color)
	reset()


##################################################################################
#
#	Fast LSA using random projection
#
##################################################################################

sn <- 2000
sp <-  300
k  <-   30

# --- structure
U  <- matrix( rnorm(sn*k), nrow=sn, ncol=k  )
D  <- 1+rchisq(k, df=2)
Vt <- matrix( rnorm(k*sp), nrow=k, ncol=sp )

X <- (U %*% (D * Vt)) + matrix( rnorm(sn*sp), nrow=sn, ncol=sp )

udv <- svd(X); u <- udv$u; plot(udv$d,log="xy")

fast.svd <- function(X, k) {
	u <- X %*% matrix(rnorm(ncol(X)*(k+5)),nrow=ncol(X),ncol=k+5)
	return (u[,1:k])
}

u.fast <- fast.svd(X,20)
plot(u[,1],u.fast[,1])


##################################################################################
#
#	 SVD via alternating least squares ALS 
#
##################################################################################


normalize <- function(x) { return (x / sqrt(x%*%x)) }

als.svd <- function(X, max.it = 5) {
	ss0 <- sum(X^2)
	u <- normalize(rowMeans(X)) # init with row means
	it <- 0;
	while(it < max.it) {
		it <- it + 1
		v <- normalize(drop(u %*%   X))   # u'u=1
		u <-           drop(v %*% t(X))   # v'v=1
		s <- drop(sqrt( u%*%u ))
		u <- u/s
		ss1 <- sum((X - s*outer(u,v))^2)
		cat(ss1,"\n")
		if (((ss0-ss1)/ss0) < .001) break
		else ss0 <- ss1
	}
	return(list(d=s,u=u,v=v))	
}


X <- matrix( rnorm(1000*200, mean=3), nrow=1000, ncol=200)

cs <- colSums(X)
rs <- rowSums(X)
X <- 1/sqrt(rs) * t( 1/sqrt(cs) * t(X) )


uv <- als.svd(X)
plot(sqrt(rs)/sqrt(sum(rs)), uv$u)
plot(sqrt(cs)/sqrt(sum(cs)), uv$v)

udv <- svd(X)
c(udv$d[1], uv$d)
plot(udv$u[,1], -uv$u, xlab="SVD", ylab="ALS", main="First U"); abline(a=0,b=1,col="gray")
plot(udv$v[,1], -uv$v, xlab="SVD", ylab="ALS", main="First V"); abline(a=0,b=1,col="gray")



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


