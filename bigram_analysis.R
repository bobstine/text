##################################################################################
# 
# Function definitions
#
##################################################################################

reset <- function() {
	# par(mfrow=c(1,1), mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)      # default
	par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,2.5,2,1)+0.1)  # bottom left top right
	}
	
# --- run the following to initialize
reset()  # sets up plot


#######################################################################################
# Spectrum of bigram matrix
#######################################################################################
spectrum <- function() {
# ------------------------------------------------------------------------------------
#   compare exact and random project singular vectors
#		exact spectrum of raw bigram is very weird (rounding?)
#      only first has high correlation, rest are occasionally correlated
#      high cancor up to last 20% of components

city <- "ChicagoOld3"
path <- paste("/Users/bob/C/text/text_src/temp/",city,"/", sep="");

ds <- scan(paste(path,"UDV_single.d.txt",sep=""))
dd <- scan(paste(path,"UDV_double.d.txt",sep=""))
par(mfrow=c(2,1))
	plot(dd, xlab="Index of Singular Value", ylab="Singular Value", log="y", 
		main="Singular Values of Raw Bigram Matrix")
	lines(ds, col="red")
	plot(dd,ds, xlab="Double Precision", ylab="Single Precision", log="xy", 
		main="Singular Values of Unweighted Bigram Matrix")
reset();

# --- log/log model for decay of singular values, spectrum of bigram matrix
all <- 1:2500
plot(dd[all], log="xy", ylab="Leading Singular Values, Bigram Matrix", xlab="Index")
y <- log(dd[all]); x <- log(all)
summary(regr <- lm(y ~ x ))          # overall fit = 12.31 - 1.55 log(i) 
lines(exp(x),exp(fitted.values(regr)), col="blue")
x2 <- (x-mean(x))^2
summary(regr <- lm(y ~ x + x2))      # 14.33 - 1.823 log(i) - 0.145 (log(i)-mean)^2
lines(exp(x),exp(fitted.values(regr)), col="red")
i <- 100:500
y <- log(dd[i]); x <- log(i)
summary(regr <- lm(y ~ x ))          # overall fit = 12.31 - 1.55 log(i) 
lines(exp(x),exp(fitted.values(regr)), col="magenta")

}

#######################################################################################
# Centroid variables from bigram matrix
#######################################################################################
centroid <- function() {

file     <- paste(path,"bigram_",nProj,"_raw.txt", sep="")
Bigram.n <- read.table(file, header=TRUE); dim(Bigram.n)
file     <- paste(path,"bigram_",nProj,"_exact.txt", sep="")
Bigram.e <- read.table(file, header=TRUE); dim(Bigram.n)

j<-3;
plot(Bigram.n[,j], Bigram.e[,j], 
	xlab=paste("Random Projection, Component",j,sep=" "), ylab="Exact Singular Vector")
abline(a=0,b=1, col="gray")
abline(h=0,col="gray", lty=3); abline(v=0,col="gray", lty=3); 

drawit <- function(k) {
	cca <- cancor(Bigram.n[,1:k], Bigram.e[,1:k])
	plot(1:k, cca$cor, xlab="Dimensions", ylab="Canonical Correlation")
	abline(v=0.8*k,col="gray")
	}

par(mfrow=c(3,1))
	drawit(100);
	drawit(200);
	drawit(400);
reset()
}


# ------------------------------------------------------------------------------------
#   Eigenwords produced by different weightings of the bigram matrix (random projection)

epath    <- "/Users/bob/C/text/text_src/temp/bigram/"
U.raw   <- as.matrix(read.table(paste(epath,"u_raw.txt", sep=""), header=TRUE)); dim(U.raw)
U.sym   <- as.matrix(read.table(paste(epath,"u_sym.txt", sep=""), header=TRUE)); dim(U.sym)
U.left  <- as.matrix(read.table(paste(epath,"u_left.txt", sep=""), header=TRUE)); dim(U.left)
U.right <- as.matrix(read.table(paste(epath,"u_right.txt", sep=""), header=TRUE)); dim(U.right)

V.raw   <- as.matrix(read.table(paste(epath,"v_raw.txt", sep=""), header=TRUE)); dim(V.raw)
V.sym   <- as.matrix(read.table(paste(epath,"v_sym.txt", sep=""), header=TRUE)); dim(V.sym)
V.left  <- as.matrix(read.table(paste(epath,"v_left.txt", sep=""), header=TRUE)); dim(V.left)
V.right <- as.matrix(read.table(paste(epath,"v_right.txt", sep=""), header=TRUE)); dim(V.right)

#     SS of the left bigram singular vectors

ss <- function(x) {  sum(x*x)  }

par(mfrow=c(2,2))
	plot(y<-apply(U.raw,2,ss), log="y", main="No Normalization",
		ylab="Sum of Squares", xlab="Left Eigenwords")
	plot(y<-apply(U.sym,2,ss), log="y", main="Symmetric Normalization",
		ylab="Sum of Squares", xlab="Left Eigenwords")
	plot(y<-apply(U.left,2,ss), log="y", main="Left Normalization",
		ylab="Sum of Squares", xlab="Left Eigenwords")
	plot(y<-apply(U.right,2,ss), log="y", main="Right Normalization",
		ylab="Sum of Squares", xlab="Left Eigenwords")
reset()


#     SS of the right bigram singular vectors
par(mfrow=c(2,2))
	plot(y<-apply(V.raw,2,ss), log="y", main="No Normalization",
		ylab="Sum of Squares", xlab="Right Eigenwords")
	plot(y<-apply(V.sym,2,ss), log="y", main="Symmetric Normalization",
		ylab="Sum of Squares", xlab="Right Eigenwords")
	plot(y<-apply(V.left,2,ss), log="y", main="Left Normalization",
		ylab="Sum of Squares", xlab="Left Eigenwords")
	plot(y<-apply(V.right,2,ss), log="y", main="Right Normalization",
		ylab="Sum of Squares", xlab="Left Eigenwords")
reset()


# ------------------------------------------------------------------------------------
#   look at variation of centroids formed by documents in eigenword space

file     <- paste("/Users/bob/C/text/text_src/temp/",city,"bigram_",nProj,"_sym.txt", sep="")
Bigram.s <- read.table(file, header=TRUE); dim(Bigram.s)
file     <- paste("/Users/bob/C/text/text_src/temp/",city,"bigram_",nProj,"_lhs.txt", sep="")
Bigram.l <- read.table(file, header=TRUE); dim(Bigram.l)
file     <- paste("/Users/bob/C/text/text_src/temp/",city,"bigram_",nProj,"_rhs.txt", sep="")
Bigram.r <- read.table(file, header=TRUE); dim(Bigram.r)

ss <- function(x) {  sum(x*x)  }

#     SS of the bigram singular vectors
par(mfrow=c(2,1))
	plot(y<-apply(as.matrix(Bigram.n[,1:1500]),2,ss), log="xy", main="No Normalization",
		ylab="Sum of Squares", xlab="Left Components, Bigram Predictor Sequence")
	plot(y<-apply(as.matrix(Bigram.n[,1501:3000]),2,ss), log="xy",
		ylab="Sum of Squares", xlab="Right Components, Bigram Predictor Sequence")
reset()

par(mfrow=c(2,1))
	plot(y<-apply(as.matrix(Bigram.s[,1:1500]),2,ss), log="xy", main="Symmetric Normalization",
		ylab="Sum of Squares", xlab="Left Components, Bigram Predictor Sequence")
	plot(y<-apply(as.matrix(Bigram.s[,1501:3000]),2,ss), log="xy",
		ylab="Sum of Squares", xlab="Right Components, Bigram Predictor Sequence")
reset()

par(mfrow=c(2,1))
	plot(y<-apply(as.matrix(Bigram.l[,1:1500]),2,ss), log="xy", main="Left Normalization (row sum 1)",
		ylab="Sum of Squares", xlab="Left Components, Bigram Predictor Sequence")
	plot(y<-apply(as.matrix(Bigram.l[,1501:3000]),2,ss), log="xy",
		ylab="Sum of Squares", xlab="Right Components, Bigram Predictor Sequence")
reset()

par(mfrow=c(2,1))
	plot(y<-apply(as.matrix(Bigram.r[,1:1500]),2,ss), log="xy", main="Right Normalization (col sum 1)",
		ylab="Sum of Squares", xlab="Left Components, Bigram Predictor Sequence")
	plot(y<-apply(as.matrix(Bigram.r[,1501:3000]),2,ss), log="xy",
		ylab="Sum of Squares", xlab="Right Components, Bigram Predictor Sequence")
reset()


abline(v=  10, col="gray"); text(  10, 0.05, round(ss(Bigram[,  10]),4))    
abline(v=1000, col="gray"); text(1000, 0.05, round(ss(Bigram[,1000]),5))

x <- 1:1500
r <- lm(I(log(y)) ~ I(log(x))); summary(r)



# --- To see how the variance of the centroids shrinks, we need to see the doc/word matrix
#     and form the centroids.  (need to remove """ in source first)
#     Need U matrix from above to compute

W <- as.matrix(read.table(paste(path,"w5708.txt",sep=""), header=TRUE)); dim(W)    # 7384 x 5708
U <- U.raw; dim(U)                                                                 # 5708 x 500

B <- W %*% U

}
