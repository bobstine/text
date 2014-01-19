

reset <- function() {
	# par(mfrow=c(1,1), mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)      # default
	par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,2.5,2,1)+0.1)  # bottom left top right
	}



##################################################################################
# linear lsa results
##################################################################################

path <- "/Users/bob/C/text/text_src/temp/ChicagoOld3/"
YM <- as.matrix(read.table( paste(path,"lsa_ym.txt",sep=""), header=TRUE, as.is=TRUE))
X  <- as.matrix(read.table( paste(path,"lsa_1500.txt",sep=""), header=TRUE, as.is=TRUE))


Y <-  YM[,1]    # log price
M <-  YM[,2]    # length of documents
lM <- log(M)
Xw <- X
for(i in 1:nrow(Xw)) { Xw[i,] <- sqrt(YM[i,2]) * X[i,]; }
#																Adjusted R2
#                                                               W'[D W1]
#														words	counts		sym wts		
sr <- summary(r <- lm(Y ~ M + lM + X [,1: 500])); sr  #	0.5809	0.6422		0.5724			
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1: 500])); srw # 		0.6509		0.5831

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1000])); sr  #	0.6311	0.6925           
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1:1000])); srw # 		0.6929

sr <- summary(r <- lm(Y ~ M + lM + X [,1:1500])); sr  #	0.6699	0.7244		0.6896 
srw<- summary(r <- lm(Y ~ M + lM + Xw[,1:1500])); srw # 		0.7186		0.6863

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

