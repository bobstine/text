
source("/Users/bob/C/text/functions.R")

# Notes:
#
#		Quadratic not so useful: maybe just too sparse?
#
##################################################################################
# linear lsa results
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

name <- "Quadratic"
X  <- as.matrix(read.table( paste(path,"lsaq_250.txt",sep=""),     header=TRUE, as.is=TRUE))


# only useful for raw counts
Xw <- X
for(i in 1:nrow(Xw)) { Xw[i,] <- X[i,] / sqrt(YM[i,2]) }

#																Adjusted R2
#														words		Raw		Sqrt	1/ni	Quad, Raw		
sr <- summary(r <- lm(Y ~ M + lM + X [,1: 250])); sr  #				0.576	0.602	0.573	0.385
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1: 250])); srw #				0.600					0.405

sr <- summary(r <- lm(Y ~ M + lM + X [,1: 500])); sr  #	0.5809		0.612	0.640	0.628
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1: 500])); srw #				0.636

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1000])); sr  #	0.6311	    0.655	0.686	0.683
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1:1000])); srw #				0.682

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1500])); sr  #	0.6699		0.682	0.714	0.715
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1:1500])); srw #				0.709


plot(coefficients(sr)[-1,3],coefficients(srw)[-1,3])
text(coefficients(sr)[ 2,3],coefficients(srw)[ 2,3], "M")


par(mfrow=c(1,2))    # regrW.pdf
	y <- abs(coefficients(sr)[-(1:3),3])
	x <- 1:length(y)
	plot(x,y, 
		xlab=name, ylab="|t|", main="")
		abline(h=-qnorm(.025/(nProj/2)), col="gray", lty=3)
		abline(h=sqrt(2/pi), col="cyan")
		lines(lowess(x,y,f=0.3), col="red")
	half.normal.plot(y,height=5)
reset()


##################################################################################
# quadratc lsa comparison
##################################################################################

# assume X is linear, Xq is quadratic


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

