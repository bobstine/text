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
# Comparison of different normalizations of the bigram matrix
#######################################################################################
bigram.analysis <- function() { 

# ------------------------------------------------------------------------------------
#   compare exact and random project singular vectors
#		exact spectrum of raw bigram is very weird (rounding?)
#      only first has high correlation, rest are occasionally correlated
#      high cancor up to last 20% of components

city <- "ChicagoOld3"

ds <- scan(paste("/Users/bob/C/text/text_src/temp/",city,"UDV_single.d.txt",sep=""))
dd <- scan(paste("/Users/bob/C/text/text_src/temp/",city,"UDV_double.d.txt",sep=""))
par(mfrow=c(2,1))
	plot(dd, xlab="Index of Singular Value", ylab="Singular Value", log="y", 
		main="Singular Values of Raw Bigram Matrix")
	lines(ds, col="red")
	plot(dd,ds, xlab="Double Precision", ylab="Single Precision", log="xy", 
		main="Singular Values of Unweighted Bigram Matrix")
reset();

file     <- paste("/Users/bob/C/text/text_src/temp/",city,"bigram_",nProj,"_raw.txt", sep="")
Bigram.n <- read.table(file, header=TRUE); dim(Bigram.n)
file     <- paste("/Users/bob/C/text/text_src/temp/",city,"bigram_",nProj,"_exact.txt", sep="")
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


# ------------------------------------------------------------------------------------
#   Eigenwords produced by different weightings of the bigram matrix (random projection)

path    <- "/Users/bob/C/text/text_src/temp/bigram/"
U.raw   <- as.matrix(read.table(paste(path,"u_raw.txt", sep=""), header=TRUE)); dim(U.raw)
U.sym   <- as.matrix(read.table(paste(path,"u_sym.txt", sep=""), header=TRUE)); dim(U.sym)
U.left  <- as.matrix(read.table(paste(path,"u_left.txt", sep=""), header=TRUE)); dim(U.left)
U.right <- as.matrix(read.table(paste(path,"u_right.txt", sep=""), header=TRUE)); dim(U.right)

V.raw   <- as.matrix(read.table(paste(path,"v_raw.txt", sep=""), header=TRUE)); dim(V.raw)
V.sym   <- as.matrix(read.table(paste(path,"v_sym.txt", sep=""), header=TRUE)); dim(V.sym)
V.left  <- as.matrix(read.table(paste(path,"v_left.txt", sep=""), header=TRUE)); dim(V.left)
V.right <- as.matrix(read.table(paste(path,"v_right.txt", sep=""), header=TRUE)); dim(V.right)


ss <- function(x) {  sum(x*x)  }

#     SS of the left bigram singular vectors
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



