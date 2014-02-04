
source("/Users/bob/C/text/functions.R")

# Notes:
#
#		Quadratic not so useful: maybe just too sparse?
#


##################################################################################
# linear and quadratic lsa results, real estate
##################################################################################

# --- response, weighted
path <- "/Users/bob/C/text/text_src/temp/ChicagoOld3/"
YM <- as.matrix(read.table( paste(path,"lsa_ym.txt",sep=""), header=TRUE, as.is=TRUE))

Y <-  YM[,1]    # log price
M <-  YM[,2]    # length of documents
lM <- log(M)

# --- versions
name <- "Linear, no wts"
X  <- as.matrix(read.table( paste(path,"lsa_raw_1500.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, sqrt"
X  <- as.matrix(read.table( paste(path,"lsa_sqrt_1500.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, recip"  # conditional prob
X  <- as.matrix(read.table( paste(path,"lsa_recip_1500.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, cca" 
X  <- as.matrix(read.table( paste(path,"lsa_cca_1500.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, tf-idf"  # conditional prob
X  <- as.matrix(read.table( paste(path,"lsa_tfidf_1500.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Quadratic"
Xq <- as.matrix(read.table( paste(path,"lsaq_raw_250_p0.txt",sep=""), header=TRUE, as.is=TRUE))


#																Adjusted R2
#														words	Raw		Sqrt	1/ni	cca		tfidf	quad		
sr <- summary(r <- lm(Y ~ M + lM + X [,1: 250])); sr  #			0.576	0.602	0.573	0.582	0.592	0.390

sr <- summary(r <- lm(Y ~ M + lM + X [,1: 500])); sr  #	0.581	0.612	0.640	0.628	0.639	0.628

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1000])); sr  #	0.631	0.655	0.686	0.683	0.697	0.664

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1500])); sr  #	0.670	0.682	0.714	0.715	0.722	0.686


plot(coefficients(sr)[-1,3],coefficients(srw)[-1,3])
text(coefficients(sr)[ 2,3],coefficients(srw)[ 2,3], "M")


quartz(width=6.5,height=3); reset()
par(mfrow=c(1,2))                                   #        regrW.pdf
	y <- abs(coefficients(sr)[-(1:7),3])
	x <- 1:length(y)
	plot(x,y, cex=0.25, col="darkgray",
		xlab=name, ylab="|t|", main="")
		abline(h=-qnorm(.025/length(y)), col="black", lty=4)
		abline(h=sqrt(2/pi), col="black", lty=2)
		smth <- loess(y~x,span=0.5)
		lines(predict(smth), col="red")
	half.normal.plot(y,height=5)
reset()


##################################################################################
# quadratic lsa comparison
##################################################################################

#  X is linear, Xq is quadratic (no weights)
#  Not much better than just 500 linear
sr <- summary(r <- lm(Y ~ M + lM + X[,1: 250] + Xq[,1:250])); sr     # 0.593

# nearly linear decay
cca <- cancor(X[,1:ncol(Xq)],Xq)
plot(cca$cor)

# weights are pretty random on quad components, decay on xcoef
plot(cca$ycoef[,1])

# add to regression												Linear alone	With Quad
sr <- summary(r <- lm(Y ~ M + lM + X [,1: 250]     )); sr  #		0.576
sr <- summary(r <- lm(Y ~ M + lM + X [,1: 250] + Xq)); sr  #					0.594


##################################################################################
# Raw word regression
##################################################################################

W <- as.matrix(read.table("/Users/bob/C/text/text_src/temp/ChicagoOld3/w5708.txt", header=TRUE)); dim(W)

sr <- summary(r <- lm(Y ~ M + lM + W [,1: 500])); sr  #	0.5809
sr <- summary(r <- lm(Y ~ M + lM + W [,1:1000])); sr  #	0.6311
sr <- summary(r <- lm(Y ~ M + lM + W [,1:1500])); sr  #	0.6699
sr <- summary(r <- lm(Y ~ M + lM + W [,1:2000])); sr  #	0.6875

# compute tf-idf (remove EOL char)
W <- W[,-7]
term.counts <- apply(W>0,2,sum)
(W[1:5,] * log(nrow(W)/term.counts))[,1:5]


##################################################################################
# linear lsa results, wine
##################################################################################

# --- response, weighted  (n=20,873)
path <- "/Users/bob/C/text/text_src/temp/wine/"
YM <- as.matrix(read.table( paste(path,"lsa_ym.txt",sep=""), header=TRUE, as.is=TRUE))

Y <-  YM[,1]    # log price
M <-  YM[,2]    # length of tasting notes
lM <- log(M)

# --- versions
name <- "Linear, no wts"
X  <- as.matrix(read.table( paste(path,"lsa_raw_1000.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, sqrt"
X  <- as.matrix(read.table( paste(path,"lsa_sqrt_1000.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, recip"  # conditional prob
X  <- as.matrix(read.table( paste(path,"lsa_recip_1000.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, CCA" 
X  <- as.matrix(read.table( paste(path,"lsa_cca_1000.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Linear, tf-idf"  # conditional prob
X  <- as.matrix(read.table( paste(path,"lsa_tfidf_1000.txt",sep=""), header=TRUE, as.is=TRUE))

name <- "Quadratic"
X  <- as.matrix(read.table( paste(path,"lsaq_250.txt",sep=""),     header=TRUE, as.is=TRUE))


#																Adjusted R2
#														Raw		Sqrt	1/ni	cca		tfidf	Quad		
sr <- summary(r <- lm(Y ~ M + lM + X [,1: 250])); sr  #	0.673	0.675	0.663	0.650	0.685

sr <- summary(r <- lm(Y ~ M + lM + X [,1: 500])); sr  #	0.706	0.708	0.698	0.689	0.710

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1000])); sr  #	0.733	0.737	0.729	0.712	0.732


quartz(width=6.5,height=3); reset()
par(mfrow=c(1,2))                                   #        regrW.pdf
	y <- abs(coefficients(sr)[-(1:7),3])
	x <- 1:length(y)
	plot(x,y, cex=0.25, col="darkgray",
		xlab=name, ylab="|t|", main="")
		abline(h=-qnorm(.025/length(y)), col="black", lty=4)
		abline(h=sqrt(2/pi), col="black", lty=2)
		smth <- loess(y~x,span=0.5)
		lines(predict(smth), col="red")
	half.normal.plot(y,height=5)
reset()
