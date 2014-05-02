
source("/Users/bob/C/text/functions.R")

# KL distance with avg truths
kl.dist<- function (p1,p2) {
	i1 <- which(p1>0); i2 <- which(p2>0)
	-0.5 * (sum(log( p2[i1]/p1[i1] ) * p1[i1]) + sum(log( p1[i2]/p2[i2] ) * p2[i2]))
	}

dither <- function(x) { x + rnorm(length(x),sd=0.05*sd(x)) }

# n<-20; kl.dist( dbinom(0:n, n, prob=0.8), dbinom(0:n, n, prob=0.3)  )
#        kl.dist(P[1,],P[2,])


##################################################################################
#
# Simulate data from topic model, Poisson version
#
#    simulated response is lognormal (eg, log price)
#
##################################################################################

	
	n.doc  <- 5000						#	number observed documents/listings  
	K      <-   10						#	number possible attributes (topics)
	
	mu <- rgamma(K, shape=2, scale=2)	#	'true' value of attributes
	# hist(mu, main="Simulated True Contributions to Log Prices", breaks=20)		
	
	Y <- rnorm(n.doc,mean=12.2,sd=0.8)	#	initial seed log price (revised below)
	
	A <-matrix(0,nrow=n.doc,ncol=K)		#	A identifies attributes in each doc, fuzzed
	b <- rbinom(n.doc,1,0.5)  			#	coin tosses
	for(i in 1:nrow(A)) {
		permute <- sample(1:K,K,replace=FALSE)
		cs <- cumsum(mu[permute])
		nk <- which(Y[i] <= cs)[1]
		if((nk > 1) & (b[i]==1)) nk <- nk-1
		A[i,permute[1:nk]]<-1; Y[i] <- cs[nk]
	}
	A.sum <- apply(A,1,sum)
	# hist(A.sum); 
	cat("Mean number of latent attributes", mean(A.sum),"\n"); fivenum(A.sum)
	
	# --- model fits perfectly if regress Y on columns of A, so add noise
	pct.noise <- 0.0
	sigma <- sqrt(var(Y)*pct.noise/(1-pct.noise))
	cat ("sigma = " , sigma,"\n")
	Y.obs <- Y + sqrt(var(Y)*pct.noise/(1-pct.noise)) * rnorm(length(Y))
	# rescale SD
	Y.obs <- mean(Y.obs) + (Y.obs-mean(Y.obs))*sd(logPrice)/sd(Y.obs)
	cat("mean and sd of Y.obs", mean(Y.obs),sd(Y.obs),"\n")
	# plot(density(logPrice), main="", xlab="log price", ylim=c(0,0.4))
	plot(density(Y.obs), col="gray")
	
	summary(r <- lm( Y.obs ~ A ))
	plot(mu,coefficients(r)[-1])
	plot(Y.obs ~ A.sum, xlab="Number Latent Attributes", ylab="Log Price"); 
	cat("Squared attribute corr: ",cor(Y,A.sum)^2,", logged ",cor(Y, log(A.sum))^2,"\n")



	n.specific	<- 1000
	n.common 	<-  200
	n.vocab <- n.common + n.specific					# size of vocabulary
	m.topic <- floor(n.specific/K)						# word types per topic
	alpha   <- 0.1  									# dirichlet parameter within topic words
	P       <- matrix(0, nrow=K, ncol=n.vocab)			# only *specific* components
	for(k in 1:K) {  
		j <- (n.common + 1+(k-1)*m.topic):(n.common + k * m.topic) 														
		P[k,j] <- rdirichlet(rep(alpha,m.topic))
	}
	apply(P,1,sum)[1:4]										# prob dist so rows sum to 1
	plot(P[1,]); points(P[2,],col="red"); points(P[4,],col="blue")
	par(mfrow=c(2,2)); 										# only Zipf has common support
		plot(P[1,],P[2,]); plot(P[3,],P[4,]); 
		plot(P[5,],P[6,]); plot(P[7,],P[8,]); 
	reset()
	
	# --- divergence distances
	c(kl.dist(P[1,],P[2,]), kl.dist(P[2,],P[3,]), kl.dist(P[3,],P[4,]), kl.dist(P[4,],P[5,]))

	# plot sparse and non sparse distributions				[ P.pdf ] 
	par(mfrow=c(1,2))
		plot(P.01[4,], P.01 [3,], log="xy",cex=0.5,xlab=expression(P[1]),ylab=expression(P[2]),
			main=expression(alpha~"=0.010"), cex.main=0.8)
		plot(P.025[4,],P.025[2,], log="xy",cex=0.5,xlab=expression(P[1]), ylab=expression(P[2]),
			main=expression(alpha~"=0.025"), cex.main=0.8)
	reset()

	# --- generate doc/word matrix W
	W 		<- matrix(0, nrow=n.doc, ncol=n.vocab)
	zipf	<- 1/(1:n.common); zipf <- zipf/sum(zipf)		# Zipf distribution
	lambda.s <- 10											# expected word tokens per topic
	one		<- rep(1,n.common)
	for(i in 1:n.doc) { 
		jj <- which(A[i,]==1)
		nt <- length(jj)									# num topics in doc i
		# n.words <- rep(lambda.s,nt);
		n.words	<- pmax(1,rpois(nt, lambda.s))
		W[i,1:n.common] <- rpois(one, nt*lambda.s*zipf)	# common words
		for (j in 1:length(jj)) {
			indx <- sample(1:n.vocab, n.words[j], replace=FALSE, prob=P[jj[j],]);
			W[i,indx] <- W[i,indx] + 1
		}
	}
	cat("Max word frequencies in first 5 docs", apply(W,1,max)[1:5])
	doc.len <- apply(W,1,sum)						# check document lengths
	c(mean(doc.len), sd(doc.len)); hist(doc.len, breaks=25)
	
	# --- cca sanity check
	j <- which(P[1,]>0.01); j
	plot(dither(rowSums(W[,j])), dither(A[,1]))

	# --- CCA analysis between attributes and words
	keep <-	which(0 < (col.sd <- apply(W,2,sd))); length(keep)
	cca <- cancor(A,W[,keep],xcenter=F,ycenter=F)
	plot(cca$cor)
	plot(cca$xcoef[,3])
	
	k <- which.max(abs(cca$xcoef[,3])); k #29
	i <- which(P[k,keep]>0.001)
	plot(P[k,keep][i], cca$ycoef[,2][i])
	max(abs(cor(P[k,keep],cca$ycoef)))
	
	
	# qqplot(doc.len,nTokens); abline(0,1,col="red")# matches but for tail [need data]

	# --- check vocab for Zipf distribution         [zipfcorr.pdf]
	# quartz(width=7.5,height=3); reset()
	par(mfrow=c(1,2))
		freq <- apply(W,2,sum)
		sort.freq <- sort(freq[freq>0],decreasing=TRUE) 
		lf <- log(sort.freq); lr <- log(1:length(sort.freq))
		plot(sort.freq,log="xy", xlab="Word Rank",ylab="Frequency")
		lf <- lf[1:500]; lr <- lr[1:500]
		regr <- lm(lf ~ lr); coefficients(summary(regr)); 
		lines(exp(predict(regr)), col="red")
		cat(sum(freq==0)," unused words\n");
	
	# --- regress log price on length
		plot(Y ~ doc.len, xlab="Document Length", ylab="Response", );  
		fit <-loess(Y ~ doc.len, span=0.4)    			# nothing nonlinear
		o <- order(doc.len)
		lines(doc.len[o], predict(fit)[o],col="red")	# some little nonlinear
		r <- lm(Y ~ doc.len); 
		cat("Sqr corr with doc length:",cor(Y,doc.len)^2," log ",cor(Y,log(doc.len))^2,"\n")
		# lines(doc.len[o], predict(r)[o], col="blue")
	reset()

#######################################################################################
#
#   Fitted regressions
#
# 		Amazing effect...
#		R2 in regr on words is 1 if fixed number per topic (know how many topics), but 
# 		R2 drops to about 0.75 if Poisson number of words per topic
#
#######################################################################################

	# --- Build ordered words (word types in order of overall frequency, dropping 0)
	freq <- apply(W,2,sum)
	o <- order(freq,decreasing=TRUE); cat("Most common words: ", o[1:20], "\n"); dim(W.ordered)
	o <- o[freq[o]>0]
	W.ordered <- W[,o]
	

	# saturated model to get max possible R2 with words + len
	sr <- summary(regr <- lm(Y.obs ~ doc.len + W.ordered));  sr
	coef.summary.plot(sr, "Word Frequencies", omit=2)

	# most common 250 words
	# sr <- summary(regr <- lm(Y.obs ~ doc.len + W.ordered[,1:500]));  sr
	# quartz(width=6.5,height=3); reset();
	# 
	
	# --- LSA analysis, raw counts
	udv <- svd(W.ordered)
	plot(udv$d[1:200], log="xy")
	U <- udv$u	
	
	# --- LSA analysis, partially scaled (cols)
	w.freq <- apply(W.ordered,2,sum)
	W.cols <- t( (1/sqrt(w.freq)) * t(W.ordered))
	udv.cols <- svd(W.cols)
	plot(udv.cols$d[1:200], log="y")
	U.cols <- udv.cols$u

	# --- LSA analysis, partially scaled (rows)
	W.rows <- (1/sqrt(doc.len)) * W.ordered
	udv.rows <- svd(W.rows)
	plot(udv.rows$d[1:200], log="y")
	U.rows <- udv.rows$u

	# --- LSA analysis, CCA scaled
	W.cca <- (1/sqrt(doc.len)) * W.ordered
	w.freq <- apply(W.ordered,2,sum)
	W.cca <- t( (1/sqrt(w.freq)) * t(W.cca))
	udv.cca <- svd(W.cca)
	plot(udv.cca$d[1:200], log="y")
	U.cca <- udv.cca$u	
		
		
	# --- Plot spectra			[ spectra.pdf, spectra2.pdf ]
	i <- 2:100; x <- log(i)
    plot(y<-udv$d[i], xlab="Component",ylab="Singular Value",log="y",type="l", ylim=c(30,105))
    y <- log(y); coefficients(lm(y ~ x))        	# -0.6     
    lines( 10 * (y<-udv.rows$d[i]), lty=4)  		# dot dash               
    y <- log(y); coefficients(lm(y ~ x))        	# -0.6     
    lines( 50 * (y<-udv.cols$d[i]), lty=3) 		# short dash
    y <- log(y); coefficients(lm(y ~ x))        	# -0.25    
    lines(500 * (y<-udv.cca $d[i]), lty=5)     		# long/short dash              
    y <- log(y); coefficients(lm(y ~ x))        	# -0.2     

	# --- show gap
	i <-25:75
	plot( i,50 * (udv.cols$d[i]), cex=0.5,pch=3,ylim=c(73,99), 	# plus
		ylab="Singular Value", xlab="Component")
    points(i,500 * (udv.cca $d[i]), pch=4, cex=0.5)     			# times
 
	# only first 10 are evident with Poisson variation in length but no zipf
	# nicely spread out over the first 50 with the addition of the Zipf
	lsa.sr <- summary(lm(Y.obs ~ doc.len + U[,1:250])); lsa.sr
	coef.summary.plot(lsa.sr, "Unscaled LSA Variables", omit=1)

	lsa.rows.sr <- summary(lm(Y.obs ~ doc.len + U.rows[,1:250])); lsa.rows.sr
	coef.summary.plot(lsa.rows.sr, "Column-scaled LSA Variables", omit=1)
	summary(lm(Y.obs ~ doc.len +  U.rows[,1:50]))

	lsa.cols.sr <- summary(lm(Y.obs ~ doc.len + U.cols[,1:250])); lsa.cols.sr
	coef.summary.plot(lsa.cols.sr, "Column-scaled LSA Variables", omit=1)
	summary(lm(Y.obs ~ doc.len +  U.cols[,1:50]))

	lsa.cca.sr <- summary(lm(Y.obs ~ doc.len + U.cca[,1:250])); lsa.cca.sr
	coef.summary.plot(lsa.cca.sr, "CCA Scaled LSA Variables", omit=1)
	summary(lm(Y.obs ~ doc.len + U.cca[,1:50]))

	quartz(width=6.5, height=3); reset()
	par(mfrow=c(1,2))
		coef.summary.plot(lsa.cols.sr, "Column Scaled Component", ylim=c(0,9.5),omit=2,show.qq=F)
		coef.summary.plot(lsa.cca.sr,  "CCA Scaled Component", omit=2, show.qq=F)
	reset()

	# --- CCA with columns of A
	# cca			<- cancor(A, U[,1:250], xcenter=F, ycenter=F)
	# cca.scaled	<- cancor(A, U.scaled[,1:250], xcenter=F, ycenter=F)

	

##################################################################################
#
#	My version of Dean's code
#
##################################################################################

# lambda <- matrix(rexp(n*k), nrow=n, ncol=k) %*% matrix(rexp(k*m), nrow=k, ncol=m)

rdrl <- function(n,a) {
    y <- rgamma(n, shape=a, scale=1)
    return(y / sum(y))
}

	n <- 100; 				# rows in doc-term 
	k <- 5;					# attributes
	m <- 2000;				# size of vocab
	my <- 1e-4				# plots

	A <- matrix(rbinom(n*k,size=1,prob=3/k),nrow=n,ncol=k)
	P <- matrix(0,nrow=k,ncol=m); for(i in 1:k) P[i,]<-rdrl(m,0.02)
	lambda<-A %*% P

	d <- svd(lambda)$d
	
	poisson<-matrix(0,nrow=n,ncol=m); 
	for(i in 1:n) poisson[i,]<- rpois(m,lambda[i,]) 
	
	normal<-matrix(0,nrow=n,ncol=m); 
	for(i in 1:n) normal[i,] <- rnorm(m, lambda[i,], sd=sqrt(lambda[i,]))

	r <- rowSums(poisson)
	c <- colSums(poisson)/sum(W.ordered)
	expect <- (r %*% t(c)) + 1/2
	cca <- poisson/sqrt(expect)
	
	plot(d, log="y", ylim=c(my,max(d)), pch=20)
	points(pmax(my,svd(poisson)$d), col="blue", pch=20, cex=0.5)
	points(pmax(my,svd(normal)$d), col="red")
	points(pmax(my,svd(cca)$d), col="green")
	


m <- matrix(rnorm(100),nrow=10,ncol=5)

d1 <- svd(m)$d
d2 <- svd(m %*% diag(c(1,1,10,10,10)))$d

plot(d1, d2)


##################################################################################
#
#	Zipf Experiment
#
#		Zipf over words produces power law (with lower power) spectrum in random words
#
##################################################################################

# --- generate zipf dist over k integers

k     <- 500
power <- 1
zpdf <- 1/ (1:k)^power; zpdf <- zpdf/sum(zpdf)
plot( zpdf, log = "xy" )

# --- sample n documents using this distribution over words

n <- 10000

W <- matrix(0,nrow=n,ncol=k)
doc.len <- 100

for(i in 1:nrow(W)) {
	df <- as.data.frame(table(sample(1:k, doc.len, replace=T, prob=zpdf)))
	W[i,df$Var1] <- df$Freq
	}
#		retain colummns with positive counts
cs <- colSums(W); cs
dim( W <- W[,cs>2] )


# --- spectrum of W

W <- t( t(W) - colMeans(W) )

plot( d <- svd(W)$d, log="xy", 	xlab="Component", ylab="Singular Values of W")

	# fit drops first and last
d <- d[2:floor(0.99*length(d))]
y <- log(d); x <- log(1:length(d))
r <- lm(y ~ x); r
lines(exp(x),exp(predict(r)), col="red")

u <- (svd(W)$u)[,1:3]

round( t(u) %*% u  ,3)

j <- 1;
sum( outer(u[,j],u[,j]) )
sum(u[,j]) ^2



##################################################################################
#
#	SVD Experiments
#
##################################################################################

n <- 6; m <- 3
x <- round( matrix(rnorm(n*m),nrow=n), 1 )

udv   <- svd(x)

udv.r <- svd(  diag(1:n) %*% x  )

udv.c <- svd(  x %*% diag(1:m)  )

udv$d
udv.r$d
udv.c$d

# multiply rows by diagonal changes range (norm for range) for larger dim
cancor(udv$u, udv.r$u, xcenter=F, ycenter=F)$cor
cancor(udv$u, udv.c$u, xcenter=F, ycenter=F)$cor

# these match up since must span R^3
cancor(udv$v, udv.r$v, xcenter=F, ycenter=F)$cor
cancor(udv$v, udv.c$v, xcenter=F, ycenter=F)$cor


# --- orthogonal matrix

n <- 6; m <- 3
x <- round( matrix(rnorm(n*m),nrow=n), 1 )

o <- qr.Q(qr(x))
round(o,3)
t(o) %*% o

udv <- svd(o)
round(udv$d,3)
# just permutes columns due to rounding
round(cbind(o,0,udv$u),3)
round(udv$v,3)

##################################################################################
#
#	Generate attributes randomly rather than from value
#
##################################################################################

	A.sum   <- 1 + rpois(rep(1,n.doc), 6)
	hist(A.sum); mean(A.sum); max(A.sum)						
	A <-matrix(0,nrow=n.doc,ncol=K)	#	A identifies attributes in each doc, fuzzed
	for(i in 1:n.doc) A[i,sample(1:K,A.sum[i])] <- rnorm(A.sum[i],mean=1,sd=0.2)
	Y <- A %*% mu	#  rnorm(n.doc,mean=0,sd=0)		#	perfect true model
	hist(Y, main="Dist of Log Price")

##################################################################################
#
#	Correlation between sum and count for several RVs; more skew, less corr
#
##################################################################################

	#	rnorm(k[i],mean=100,sd=10)
	
	s <-rep(0,300)
	
	k <- rpois(length(s),3)
	for(i in 1:length(s)) s[i] <- sum(rgamma(k[i],shape=0.01,scale=5))
	
	plot(k,s); cor(k,s)
	
	hist(s)

##################################################################################
#
# Simulate data from topic model, original version
#
##################################################################################
original.topic.model <- function() { 

# --- functions

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

