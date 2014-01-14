

reset <- function() {
	# par(mfrow=c(1,1), mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)      # default
	par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,2.5,2,1)+0.1)  # bottom left top right
	}


##################################################################################
# Simulate data
##################################################################################

vocab.size <- 50
n.tokens <- 1000
n.docs   <- 20
doc.len  <- n.tokens/n.docs

tokens <- sample(1:vocab.size, n.tokens, replace=TRUE)
docs   <- 1+ floor( (0:(n.tokens-1))/doc.len)


W <- matrix(0, nrow=n.tokens, ncol=vocab.size)
for(i in 1:nrow(W)) W[i,tokens[i]] <- 1

W0 <- W
W1 <- matrix(0, nrow=n.tokens, ncol=vocab.size)
W1[1:n.tokens-1,] <- W0[2:n.tokens,]


D <- matrix(0, nrow=n.tokens, ncol=n.docs)
for(i in 1:nrow(D)) D[i,docs[i]] <- 1


##################################################################################
# Compare CCA to eigenvectors and SVD
##################################################################################

cca <- cancor (W0, W1, xcenter=FALSE, ycenter=FALSE); dim(cca$xcoef)

	B01 <- t(W0) %*% W1
	B10 <- t(B01)

	F0  <- apply(W0, 2, sum)     
	F0i <- 1/F0
	F1  <- apply(W1, 2, sum)                
	F1i <- 1/F1

ev  <- eigen( diag(F0i) %*% B01 %*% diag(F1i) %*% B10 )

	H0i <- 1/sqrt(F0)
	H1i <- 1/sqrt(F1)

udv <- svd(diag(H0i) %*% B01 %*% diag(H1i) )


# --- evalues correspond to squared can cor; evectors perfectly corr
cca$cor[1:5]^2
ev$values[1:5]

j <- 2; plot(cca$xcoef[,j], ev$vectors[,j]); cor(cca$xcoef[,j], ev$vectors[,j])


# --- singular values match can cor; rescaled singular vectors perfect corr
cca$cor[1:5]
udv$d[1:5]

H <- diag(H0i)
j <- 2; plot(cca$xcoef[,j], H %*% udv$u[,j]); cor(cca$xcoef[,j], H %*% udv$u[,j])



##################################################################################
# Centroids of LSA = left evectors of SVD
##################################################################################

# lsa formed from left singular vectors of D'W, or right from W'D
udv <- svd(t(W) %*% D)

lsa <- udv$v

# alternatively form centroids of documents from left singular vectors (mean shifts coefs)
u      <- udv$u
center <- matrix(0, n.docs, ncol(u))
for(i in 1:n.docs) center[i,] <- apply(u[tokens[which(docs==i)],],2,sum)

# similarity: if docs don't have the same length, then you need to put that scaling into
# D matrix... otherwise you lose the linearity.

par(mfrow=c(2,2))
	j <- 1; plot(x<-lsa[,j], y<-center[,j]); abline(r<-lm(y~x),col="red"); round(coefficients(r)[2],2)
	j <- 2; plot(x<-lsa[,j], y<-center[,j]); abline(r<-lm(y~x),col="red"); round(coefficients(r)[2],2)
	j <- 3; plot(x<-lsa[,j], y<-center[,j]); abline(r<-lm(y~x),col="red"); round(coefficients(r)[2],2)
	j <- 4; plot(x<-lsa[,j], y<-center[,j]); abline(r<-lm(y~x),col="red"); round(coefficients(r)[2],2)
reset()

b <- rep(0,length(udv$d))
for(j in 1:length(b)) {x<-lsa[,j]; y<-center[,j]; b[j] <- coefficients(lm(y~x))[2]}
plot(b,udv$d); abline(r<-lm(udv$d ~ b), col="red"); coefficients(r)

cca <- cancor(lsa,center)

cca$cor
round(cca$xcoef,2)



##################################################################################
# Unified SVD results
##################################################################################

path <- "/Users/bob/C/text/text_src/temp/ChicagoOld3/"
YM <- as.matrix(read.table( paste(path,"unified_y_m.txt",sep=""), header=TRUE, as.is=TRUE))
X  <- as.matrix(read.table( paste(path,"unified_svd_1500.txt",sep=""), header=TRUE, as.is=TRUE))


Y <-  YM[,1]    # log price
M <-  YM[,2]    # length of documents
lM <- log(M)
Xw <- X
for(i in 1:nrow(Xw)) { Xw[i,] <- sqrt(YM[i,2]) * X[i,]; }
#																Adjusted R2
#														words	counts			sym wts
sr <- summary(r <- lm(Y ~ M + lM + X [,1: 500])); sr  #	0.5809	0.6422			0.5724
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1: 500])); srw # 		0.6509			0.5831

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1000])); sr  #	0.6311	0.6925           
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1:1000])); srw # 		0.6929

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1500])); sr  #	0.6699	0.7244			0.6896 
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1:1500])); srw # 		0.7186			0.6863

plot(coefficients(sr)[-1,3],coefficients(srw)[-1,3])
text(coefficients(sr)[ 2,3],coefficients(srw)[ 2,3], "M")


par(mfrow=c(1,2))    # regrW.pdf
	y <- abs(coefficients(sr)[-(1:2),3])
	x <- 1:length(y)
	plot(x,y, 
		xlab="Unified Predictor", ylab="|t|", main="")
		abline(h=-qnorm(.025/(nProj/2)), col="gray", lty=3)
		abline(h=sqrt(2/pi), col="cyan")
		lines(lowess(x,y,f=0.3), col="red")
	half.normal.plot(y,height=5)
reset()


##################################################################################
# Raw word regression
##################################################################################

W <- as.matrix(read.table("/Users/bob/C/text/text_src/temp/ChicagoOld3/w5708.txt", header=TRUE)); dim(W)

sr <- summary(r <- lm(Y ~ M + lM + W [,1: 500])); sr  #	0.5809
sr <- summary(r <- lm(Y ~ M + lM + W [,1:1000])); sr  #	0.6311
sr <- summary(r <- lm(Y ~ M + lM + W [,1:1500])); sr  #	0.6699
sr <- summary(r <- lm(Y ~ M + lM + W [,1:2000])); sr  #	0.6875

