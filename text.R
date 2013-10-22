##################################################################################
# 
# Function definitions
#
##################################################################################

library(xtable)    # latex tables
library(MASS)

reset <- function() {
	# par(mfrow=c(1,1), mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)      # default
	par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,2.5,2,1)+0.1)  # bottom left top right
	}

half.normal.plot <- function(y, height=2) {
	k <- length(y)          
	plot(x <- qnorm(.5+(0:(k-1))/(2*k)), y <- sort(y), 
		xlab="Normal Quantile", ylab="Sorted |t|"); 
	k <- floor(0.2*k)
	cat("Using lower 20% of cases (",k,")\n"); 
	regr <- lm(y[1:k] ~ x[1:k])
	abline(0,1, col="gray")
	abline(regr,col="red")
	text(0.25,height,paste("b =",round(coefficients(regr)[2],1)), cex=0.7)
	summary(regr)
	}
	
jitter <- function(x) { x + 0.05 * sd(x) * rnorm(length(x)) }

cross.validate.mse <- function(regr, n.folds=10,seed=23479,n.reps=1) { 	
	yx <- data.frame(y=regr$y, regr$x[,-1]) # remove intercept
	n <- nrow(yx);
	i <- rep(1:n.folds,ceiling(n/n.folds))
	if (length(i) != n) cat("Note: n is not multiple of # folds.\n");
	i <- i[1:n]
	set.seed(seed)
	mse <- rep(0,n.folds*n.reps)
	for(kk in 1:n.reps) {
		ii <- sample(i, n)   # permute fold indices
		for(fold in 1:n.folds) {
			cat(fold," ");
			train <- which(fold != ii)
			r <- lm(y ~ ., data=yx[train,])
			test  <- which(fold == ii)
			err <- yx$y[test] - predict(r, newdata=yx[test,]);
			mse[fold+(kk-1)*n.folds] <- sum(err^2)/length(test)
		}}
	mse
	}
	
show.cv <- function(regr, mse=NULL, reps=1, seed=2382) {
	if(is.null(mse)) { mse <- cross.validate.mse(regr, n.reps=reps, seed=seed)}
	hist(sqrt(mse), main=paste("Cross Validation", deparse(formula(regr))))
	red  <- (sr<-summary(regr))$sigma       ; abline(v=red , col="red", lty=3)
	red2 <- red * sqrt(1+sr$df[1]/sr$df[2]) ; abline(v=red2, col="red")
	blue <- sqrt(mean(mse))                 ; abline(v=blue, col="blue")
	n <- length(mse)
	ci <- sqrt(mean(mse)+c(-2,2)*sd(mse)/sqrt(n))
	rect(ci[1],0,ci[2],0.25,col="lightblue",border=NA)
	list(fit=red, cv=blue, mse=mse)
	}
	
# --- run the following to initialize
reset()  # sets up plot


##################################################################################
#  type counts, zipf
##################################################################################
zipf.plot <- function() { 

# --- look at type frequencies, zipf plot   (zipf.pdf)
type.cts <- sort(scan("/Users/bob/C/text/text_src/temp/type_freq.txt"), decreasing=TRUE)

x<-1:length(type.cts); y<-type.cts
zipf.data <- data.frame(list(x=x,y=y,lx=log(x),ly=log(y)))

plot(y~x, xlab="rank", ylab="frequency", log="xy", data=zipf.data)
common.words <- c(".", ",", "and", "-", "in")
text(0.9*x[1:5],0.7*y[1:5],common.words,cex=c(1,1,0.5,1,0.5))

regr<-lm(ly~lx, data=zipf.data[1:500,]); coefficients(regr)
lx <- log(x<-c(1,5000)); y <- exp(predict(regr, data.frame(lx=lx)))
lines(x,y,col="red")

}

##################################################################################
# Analysis of text regressors for real estate
##################################################################################
import.data <- function() { 

# --- analysis of regression models
nProj <- 1000

city  <- "ChicagoOld3"
file  <- paste("/Users/bob/C/text/text_src/temp/",city,"_bigram_regr.txt",sep="")

Data <- read.table(file, header=TRUE); dim(Data)

n <- nrow(Data)
price     <- Data[,"Y"]
logPrice  <- as.numeric(log(Data[,"Y"]))
nTokens   <- Data[,"m"]
logTokens <- log(nTokens)


# --- lengths (m)
mean(nTokens); fivenum(nTokens); quantile(nTokens,0.87)
boxplot(nTokens, horizontal=TRUE, xlab="Lengths of Descriptions")   # boxplot.pdf
hist(log10(nTokens))


# --- analysis of prices (thousands of $)
par(mfrow=c(1,2))                                             # prices.pdf
	y <- price
	hist(log10(price), breaks=30, main=" ",xlab="log10(Price)")
	qqnorm(log10(price), ylab="log10(Price)"); abline(a=mean(log10(y)),b=sd(log10(y)))
reset()

# --- simple models for log of prices has discontinuity 
plot(logPrice ~ logTokens) 
lines(lowess(logTokens, logPrice, f=.3), col="red")
}

##################################################################################
# Parsed variables
##################################################################################
parsed.analysis <- function() {
	
sqft  <- Data[,"SqFt"];      sqft.obs <- 0<Data[,"SqFt"]         
sqft[!sqft.obs] <- mean( sqft[sqft.obs] )
baths <- Data[,"Bathrooms"]; bath.obs <- 0<Data[,"Bathrooms"]
baths[!bath.obs] <- mean( baths[bath.obs] )
beds  <- Data[,"Bedrooms"];  beds.obs <- 0<Data[,"Bedrooms"]
beds[!beds.obs] <- mean( beds[beds.obs] )

# --- percentages missing in parsed data
(n-sum(sqft.obs))/n
(n-sum(bath.obs))/n
(n-sum(beds.obs))/n

# --- plots of the parsed explanatory variables and response           parsed.pdf
par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2,1,0))
	plot(logPrice ~ logTokens, ylab= "Log Price",  xlab="Log Number of Tokens")
	  text(1.5,6, paste("r =",round(cor(logPrice,logTokens),2)))
	  # lines(lowess(logTokens, logPrice,f=0.25), col="red")
	plot(logPrice ~ baths, ylab= "Log Price",
	  xlab="Number Bathrooms   (74% missing)", col="gray") 
	  text(5,6, paste("r =",round(cor(logPrice,baths),2)))
	points(baths[which(1==bath.obs)], logPrice[which(1==bath.obs)])  # overplot gray
	plot(logPrice ~ beds , ylab= "Log Price", 
	  xlab="Number Bedrooms  (58% missing)"  , col=c("gray","black")[1+beds.obs])   
	  text(7,6, paste("r =",round(cor(logPrice,beds),2)))
	plot(logPrice ~ I(log(sqft)),  ylab= "Log Price", 
	  xlab="Log(Sq Ft)  (94% missing)", col="gray")
	  text(2,6, paste("r =",round(cor(logPrice,log(sqft)),2)))
	points(log(sqft)[which(1==sqft.obs)], logPrice[which(1==sqft.obs)])
par(mfrow=c(1,1))

#     corr with sqft and logprice for not missing
cor( log(sqft)[sqft.obs], logPrice[sqft.obs] )
cor( baths[bath.obs], logPrice[bath.obs])

#     anomalies:  25 @ 265
table(nTokens)
ii <- which(nTokens==265); length(ii)  

# --- parsed regression fit

x.parsed. <- cbind(log(Data[,"m"]), log(sqft), sqft.obs, beds, beds.obs, baths, bath.obs)
colnames(x.parsed.)<-c("Log.m","Log.SqFt","SqFt.obs",
               "Bedrooms","Bedroom.obs","Bathrooms","Bathroom.obs")

summary(regr.parsed        <- lm(logPrice ~ x.parsed., x=TRUE, y=TRUE))

xtable(regr.parsed)

mse <- show.cv(regr.parsed,5)

}

##################################################################################
#  Raw word regression
##################################################################################
word.regression <- function () {
	

# --- plot of cum R2 statistic
r2.words <- read.table("/Users/bob/C/text/text_src/temp/w_regr_r2.txt", header=TRUE, as.is=TRUE)
r2.words[c(1,2,3,4,5,10,100,nrow(r2.words)),]

nr <- nrow(r2.words);
m <- r2.words[nr,"r2"]

plot(r2.words[,"r2"], type="l", xlab="Word Index", ylab="Cumulative R2")

plot(r2.words[,"r2"], type="l", log="x", xlab="Word Index", ylab="Cumulative R2")
abline(h=m, col="gray")
text(nr/4, m, round(m,2), col="gray")
lines( seq(1,nr,length.out=1000), seq(1/nr,nr/n,length.out=1000)  , col="red")  # null rate


d <- diff(c(0,r2.words[,"r2"]))
plot(d, xlab="Word Index", ylab="R2")

indx<-order(-d)
cbind(indx[1:10], round(d[indx][1:10],4), r2.words[indx[1:10],1])


# --- regression with W count matrix  (have to remove """ from names line)
W <- as.matrix(read.table("/Users/bob/C/text/text_src/temp/w2000.txt", header=TRUE)); dim(W)

summary( lm(logPrice ~ W[,1:100]) ); colnames(W)[1:10]

sr <- summary(  wregr <- lm(logPrice ~ W)   )
ts <- abs(coefficients(sr)[,3])
xtable(sr$coefficients[(order(-ts)[1:15]),], digits=c(0,4,4,2,4))

#  Residual standard error: 0.6818 on 5405 degrees of freedom
#  Multiple R-squared: 0.7662,	Adjusted R-squared: 0.6807 
#  F-statistic: 8.957 on 1978 and 5405 DF,  p-value: < 2.2e-16 
	
par(mfrow=c(1,2))                     # tstatRegrInd.pdf
	y <- ts[-1]             # drop intercept
	x <- 1:length(y)        # some may be singular
	threshold <- -qnorm(.025/length(y))
	plot(x,y,	xlab="Word Column in W", ylab="|t|", main="")
	abline(h=threshold, col="gray", lty=4)
	abline(h=sqrt(2/pi), col="cyan")
	lines(lowess(x,y,f=0.3), col="red")
	half.normal.plot(y)
reset()

#     number bigger than Bonferroni; indexing is a mess due to dropped columns
vars <- sapply(names( ts[ts>threshold] )[-1], function(s) substring(s,2))
colnames(x <- W[,vars])

#     regr on just those exceeding bonferroni
summary(  bregr <- lm(logPrice ~ x )  )

# --- residuals show only a vague hint of heteroscedasticity
r <- residuals(wregr); f <- fitted.values(wregr)
plot(r,f)

plot(wregr)
}

##################################################################################
# SVD variables, W
##################################################################################
lsa.analysis <- function() {
	
# --- LSA analysis from matrix W
kw <- 500
x.lsa.    <- as.matrix(Data[,paste("D",0:(kw-1), sep="")])

# write.csv(cbind(logPrice,x.lsa.), "~/Desktop/regr.csv")

m <- sqrt(Data$m)

p <- 500
sr <- summary(regr.lsa <- lm(logPrice ~ x.lsa. , x=TRUE, y=TRUE)); sr  # weights=m,

xtable(regr.lsa)

par(mfrow=c(1,2))    # regrW.pdf
	plot(	x <- 1:(nProj/2), y <- abs(coefficients(sr)[-1,3]), 
		xlab="Singular Vector of W", ylab="|t|", main="")
		abline(h=-qnorm(.025/(nProj/2)), col="gray", lty=3)
		abline(h=sqrt(2/pi), col="cyan")
		lines(lowess(x,y,f=0.3), col="red")
	half.normal.plot(y,height=5)
reset()

# --- residuals only hint at heteroscedasticity
plot(regr.lsa)
sres <- stdres(regr.lsa)
plot(nTokens, abs(sres))
lines(lowess(nTokens, abs(sres)), col="red")

# --- cross-validation
mse <- show.cv(regr.lsa, reps=20, seed=33213)  # use to compute first time

show.cv(regr.lsa, mse=save.mse$mse, reps=20, seed=33213)  # use if already computed

#     preliminary interactions
frame <- data.frame(logPrice,x.lsa.[,1:20])
br  <- lm(logPrice ~ .      , data = frame); summary(br)
br2 <- lm(logPrice ~ . + .*., data = frame); summary(br2)
anova(br,br2)
cor(fitted.values(regr.lsa), f <- fitted.values(br2))

}

##################################################################################
# SVD variables, B
##################################################################################

kb <- 500                                            # left and right
x.bigram. <- as.matrix(cbind(	Data[,paste("BL",0:(kb-1), sep="")],
									Data[,paste("BR",0:(kb-1), sep="")]  ))
summary(regr.bigram       <- lm(logPrice ~ x.bigram., x=TRUE,y=TRUE))  

par(mfrow=c(1,2))    # regrB.pdf
	y <- abs(coefficients(summary(regr.bigram))[-1,3])
	x <- rep(1:(nProj/2),2)          
	plot(x,y,xlab="Correlation Variable from Bigram", ylab="|t|", main="")
		abline(h=-qnorm(.025/(nProj/2)), col="gray", lty=3)
		abline(h=sqrt(2/pi), col="cyan")
		lines(lowess(x,y, f=0.3), col="red")
	half.normal.plot(y)
reset()

#     using just half
X <- as.matrix(Data[,paste("BL",0:(kb-1), sep="")])
summary(regr.bigram.left       <- lm(logPrice ~ X, x=TRUE,y=TRUE))  # just left

par(mfrow=c(1,2))    # regrBleft.pdf
		y <- abs(coefficients(summary(regr.bigram.left))[-1,3])
		x <- 1:length(y)
	plot(x, y, xlab="Correlation Variable from Bigram", ylab="|t|", main="")
		abline(h=-qnorm(.025/length(y)), col="gray", lty=3)
		lines(lowess(x,y), col="red")
		abline(h=sqrt(2/pi), col="cyan")
	half.normal.plot(y)
reset()


# --- CCA of bigram left/right decomp
left  <-   1        : (nProj/2)
right <- (nProj/2+1):  nProj
ccw <- cancor(x.bigram.[,left], x.bigram.[,right])
plot(ccw$cor, xlab="Index of Vector", ylab="Canonical Correlation")          # cca.pdf

cx <- x.bigram.[, left] %*% ccw$xcoef		# canonical vars
cy <- x.bigram.[,right] %*% ccw$ycoef

cor(cx[,1],cy[,1])

summary(cxregr <- lm(logPrice ~ cx )); d <- ncol(cx)

par(mfrow=c(1,2))                                              # regrBcca.pdf
	plot( x <- rep(1:d),                      
		y <- abs(coefficients(summary(cxregr))[-1,3]), 
		xlab="Canonical Variable from Bigram", ylab="|t|", main="")
		abline(h=-qnorm(.025/(nProj/2)), col="gray", lty=3)
		lines(lowess(x,y), col="red")
		abline(h=sqrt(2/pi), col="cyan")
	half.normal.plot(y, height=5)
reset()

# --- CCA of bigram left with LSA
left  <-   1        : (nProj/2)
ccw <- cancor(x.bigram.[,left], x.lsa.)
plot(ccw$cor, xlab="Index of Vector", ylab="Canonical Correlation")



##################################################################################
# Combined SVD variables
##################################################################################

#     sweep lsa from bigram variables
bigram.left  <- x.bigram.[,1:kb]
bigram.right <- x.bigram.[,(kb+1):(2*kb)]

ccw <- cancor(bigram.left, bigram.right)
cl <- bigram.left  %*% ccw$xcoef		# canonical vars
colnames(cl) <- paste("Left",1:ncol(cl), sep="_")

summary(r1 <- lm(logPrice ~ x.lsa.))                          # adj.r2 = 0.612  ChicagoOld3
summary(r2 <- lm(logPrice ~ x.lsa. + cl ))                    #          0.68
summary(r3 <- lm(logPrice ~ x.lsa. + cl + bigram.right))      #          0.70
summary(r4 <- lm(logPrice ~ x.lsa. + cl + x.bigram.right + x.parsed.)) # 0.70


y <- abs(coefficients(summary(r3))[-1,3])
d <- length(y)
par(mfrow=c(1,2))               
	plot( 
		x <- rep(1:d), y,                     
		xlab="Variable Index", ylab="|t|", main="")
		abline(h=-qnorm(.025/d), col="gray", lty=3)
		lines(lowess(x,y,f=0.3), col="red")
		abline(h=sqrt(2/pi), col="cyan")
	half.normal.plot(y, height=5)
reset()

anova(r3,r4)

r <- residuals(r4)   #   s = 0.654
f <- fitted.values(r4)

par(mfrow=c(1,2))
	plot(f,logPrice, xlab="Fitted Values", ylab="Log Price", main="Calibration Plot"); 
	abline(a=0,b=1,col="red")
	lines(lowess(f,logPrice,f=0.2),col="cyan")

	s <- 0.654
	plot(f, abs(r), xlab="Fitted Values", ylab="Absolute Residual", main="Residual Plot"); 
	abline(h=s * sqrt(2/pi), col="red")
	lines(lowess(f,abs(r),f=0.2),col="cyan")
reset()


##################################################################################
# Write SVD variables to C++
##################################################################################

write.data("/Users/bob/C/auctions/data/text/text_data.txt", 500)

write.data <- function(filename, n.cols) {
	write(n, filename)
	write.vec(logPrice, "Log Price", "role y", filename)
	write.mat(x.lsa.[,1:n.cols], "LSA"  , filename)
	write.mat( cl   [,1:n.cols], "Left" , filename)
	write.mat( cr   [,1:n.cols], "Right", filename)
	}

write.vec <- function(vector, name, desc, filename)  {
	write(name,   filename, append=TRUE)
	write(desc,   filename, append=TRUE)
	write(vector, filename, append=TRUE,ncolumns=length(vector))
	}

write.mat <- function(mat, name, filename) {
	names <- colnames(mat)
	desc  <- paste("role x stream", name)
	for (j in 1:ncol(mat)) write.vec(mat[,j], paste(name,names[j],sep="_") ,desc, filename)
	}

##################################################################################

# --- regr using CCA of two sets
summary(regr <- lm(logPrice ~ x.bigram.[,left] + x.lsa.)); 

#     these are same fits 
summary(regr  <- lm(logPrice ~ x.bigram.[,left]));
summary(regr1 <- lm(logPrice ~ cx              )); 

par(mfrow=c(1,2))                                              # regrBcca2.pdf
	y <- abs(coefficients(summary(regr))[-1,3]) 
	half.normal.plot(y, height=2)
	y1 <- abs(coefficients(summary(regr1))[-1,3])      
	half.normal.plot(y1, height=5)
reset()
plot(coefficients(regr1), coefficients(regr))

#    sweep the bigram left cca from LSA to get more orthogonal (or vice versa)
r   <- lm(x.lsa. ~ cx)
res <- residuals(r)
s   <- summary(lm(logPrice ~ cx))                      # adj.r2 = 0.63  ChicagoNew
s   <- summary(lm(logPrice ~ cx + res)); s             #          0.70
s   <- summary(lm(logPrice ~ cx + res + x.parsed.)); s #          0.705

par(mfrow=c(1,2))                                             
	d <- 2 * ncol(cx)
	plot( x <- rep(1:d),                      
		y <- abs(coefficients(s)[-1,3]), log="",
		xlab="Left CCA B, then Res LSA", ylab="|t|", main="")
		abline(h=-qnorm(.025/(nProj/2)), col="gray", lty=3)
		lines(lowess(x,y), col="red")
	half.normal.plot(y, height=5)
reset()


##################################################################################
# Comparison of regression models
##################################################################################

summary(regr.parsed        <- lm(logPrice ~ x.parsed.         ))

summary(regr.parsed.lsa    <- lm(logPrice ~ x.parsed. + x.lsa.))

summary(regr.parsed.bigram <- lm(logPrice ~ x.parsed. + x.bigram.))

summary(regr.all           <- lm(logPrice ~ x.parsed. + x.lsa. + x.bigram.))
summary(regr.all)$adj.r.squared

summary(regr.text          <- lm(logPrice ~             x.lsa. + x.bigram.))
summary(regr.text)$adj.r.squared


# --- compare nested models
anova(regr.text, regr.all)

# adding lsa/bigram to parsed + bigram/lsa
anova(regr.parsed.lsa   , regr.all)
anova(regr.parsed.bigram, regr.all)

# adding lsa/bigram to parsed
anova(regr.parsed, regr.parsed.lsa)
anova(regr.parsed, regr.parsed.bigram)

# adding parsed to bigram/lsa
anova(regr.bigram, regr.parsed.bigram)
anova(regr.lsa   , regr.parsed.lsa)


##################################################################################
# Lighthouse variables
##################################################################################

# --- estimated parsed variables in place of originals
s <- summary(r <- lm(x.parsed. ~ x.lsa. ))
for(j in 1:length(s)) {
	cat(names(s)[j]," ", s[[j]]$r.squared,"\n")  }

x.parsed.hat. <- fitted.values(r)
colnames(x.parsed.hat.) <- paste("Est",colnames(x.parsed.), sep=".")

#  pick a variable
yx <- as.data.frame(cbind(logPrice, x.parsed.))         # missing 
summary( lm(logPrice ~ Bathrooms, data=yx ) )

use <- 1 == x.parsed.[,"Bathroom.obs"]                      # few observed
yx <- as.data.frame(cbind(logPrice, x.parsed.)[use,])   # just obs
summary( lm(logPrice ~ Bathrooms, data=yx ) )

yx <- as.data.frame(cbind(logPrice, x.parsed.hat.))
plot(logPrice ~ Est.Bathrooms, data=yx )
summary( regr <- lm(logPrice ~ Est.Bathrooms, data=yx ) )
abline(regr, col="red")
text(3.5,6, paste("r=",round(sqrt(summary(regr)$r.squared),2)))

#     m is worse when estimated, but all of the others improve
j <- 2;
summary(lm(logPrice ~ x.parsed.[,j] + x.parsed.hat.[,j]))

#     all of them
summary(regr.parsed)
summary(regr.parsed.hat   <- lm(logPrice ~ x.parsed.hat.))

	
# --- canonical correlations are quite large, drop slowly
ccb <- cancor(x.lsa., x.bigram.)
plot(ccb$cor)

# --- cca of the left/right bigram variables; inverted hockey stick
#     Related to how many to keep?  
#     Clearly useful to trim left/right sets down:
#        only keep 'right' subspace that's not redundant with left
left  <-   1        : (nProj/2)
right <- (nProj/2+1):  nProj
ccw <- cancor(x.bigram.[,left], x.bigram.[,right])
plot(ccw$cor)


##################################################################################
#
# Simulate data from topic model
#
##################################################################################
                                                      topic.model <- function() { }

# --- functions

rdirichlet <- function(a) {
    y <- rgamma(length(a), a, 1)
    return(y / sum(y))
}

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


# --- Analysis of C++ results
#
#      RUN C++ here ( make dosim )
#

# --- look at type frequencies, zipf plot
type.cts <- sort(scan("/Users/bob/C/text/text_src/temp/type_freq.txt"), decreasing=TRUE)

x <- log(1:length(type.cts)); y<-log(type.cts)
plot(x,y, xlab="log rank", ylab="log frequency")        # simzipfa.pdf simzipfb.pdf
abline(regr<-lm(y[1:500]~x[1:500]),col="red"); coefficients(regr)


# --- get regression data
file  <- "/Users/bob/C/text/text_src/temp/sim_regr.txt"

Data <- read.table(file, header=TRUE); dim(Data)

nTokens  <- as.numeric(Data[,"n"])		# matches m above
resp <- Data[,"Y"]						# matches Y above

nProj <- 50

x.lsa.    <- as.matrix(Data[,paste( "D",0:(nProj/2-1), sep="")])
x.bigram. <- as.matrix(cbind(Data[,paste("BR",0:(nProj/2-1), sep="")],
							Data[,paste("BL",0:(nProj/2-1), sep="")]))

summary(regr.lsa           <- lm(Y ~ x.lsa.   ))
summary(regr.bigram        <- lm(Y ~ x.bigram.))        #   better fit

plot(fitted.values(regr.lsa),fitted.values(regr.bigram))
lines(lowess(fitted.values(regr.lsa),fitted.values(regr.bigram), f=1/3),  col="red")

cor(fitted.values(regr.lsa),fitted.values(regr.bigram)) # high corr if good fits  simfits.pdf

# --- cca shows several that are large (depending on separation)
ccb <- cancor(x.lsa., x.bigram.); plot(ccb$cor)

ccx <- x.lsa. 	%*% ccb$xcoef			# canonical vars
ccy <- x.bigram.	%*% ccb$ycoef

j <- 1; cor(ccx[,j],ccy[,j]); x <- ccx[,j]

summary( lm(Y ~ ccx) )

# --- cca of the left/right bigram variables; inverted hockey stick
#     these are the same subspaces; better count of number traits
left  <-   1        : (nProj/2)
right <- (nProj/2+1):  nProj
ccw <- cancor(x.bigram.[,left], x.bigram.[,right]); 
plot(ccw$cor, xlab="Canonical Variable", ylab="Canonical Correlation")           # simccab.pdf

cx <- x.bigram.[, left] %*% ccw$xcoef		# canonical vars
cy <- x.bigram.[,right] %*% ccw$ycoef

cor(cx[,1],cy[,1])
summary( lm(Y ~ cx) )

plot( cancor(cx,ccx)$cor )

# --- canonical corr of recovery with underlying structure
cc <- cancor(X,x.bigram.); plot(cc$cor)
cc <- cancor(X,x.lsa.); plot(cc$cor)


# --- linear separation of distributions
udv <- svd(P); plot(udv$d)


##################################################################################
# Marginal distributions of types and POS
##################################################################################

data <- readLines("/Users/bob/C/text/results/margins.txt")


types <- as.factor(scan(textConnection(data[2]), what=numeric(0), sep=","))
tabulate(types)
table(types)

pos <- scan(textConnection(data[4]), what=numeric(0), sep=" ")
pos <- pos[!is.na(pos)]

hist(log(types))

hist(log(pos))

# make sure sorted in descending order
types <- sort(types, decreasing=TRUE)

# percentage in largest types
n <- sum(types)

sapply(c(50,100,250,500,750,1000), function(k) c(k,sum(types[1:k])/n))



##################################################################################
# Comparison of cluster analysis of random projection matrix
##################################################################################

proj <- read.table("/Users/bob/Desktop/kmeans_data.txt"); dim(proj)

norm <- function(x) { sqrt(sum(x*x)) }

norm(proj[  1,1:50]);  norm(proj[  1,51:100])
norm(proj[ 10,1:50]);  norm(proj[ 10,51:100])
norm(proj[100,1:50]); norm(proj[100,51:100])

km <- kmeans(proj,200,iter.max=30)



##################################################################################
# Summary t stat plots for pure noise
##################################################################################
											noise.model <- function() {}
z <- rnorm(500)

par(mfrow=c(1,2))    # noise plots
	y <- abs(z)  
	h <- -qnorm(.025/length(y) )                           
	plot( 
		x <- rep(1:length(y)), y, 
		xlab="Random Noise", ylab="|t|", main="", ylim=c(0,1.1*h))
		abline(h=h, col="gray", lty=4)
		abline(h=sqrt(2/pi), col="cyan")
		lines(lowess(x,y), col="red")
	half.normal.plot(y, height=5)
reset()



##################################################################################
# Analysis of text regressors for wine
##################################################################################
											wine.model <- function() {}
nProj <- 500

file  <- paste("/Users/bob/C/text/text_src/temp/wine_regr.txt",sep="")

Data <- read.table(file, header=TRUE); dim(Data)

n <- nrow(Data)
rating    <- Data[,"Y"]
nTokens   <- Data[,"m"]

hist(rating)  ; mean(rating)   # ≈87, more bell-shaped
hist(nTokens) ; mean(nTokens)  # ≈42

# --- regression data
kw <- 250
x.lsa.    <- as.matrix(Data[,paste("D",0:(kw-1), sep="")]); 
dim(x.lsa.)

kb <- 250
x.bigram. <- as.matrix(cbind(	Data[,paste("BL",0:(kb-1), sep="")],
									Data[,paste("BR",0:(kb-1), sep="")]  ))
dim(x.bigram.)

#    LSA regression
regr <- lm(rating ~ x.lsa.)
s <- summary(regr); s

par(mfrow=c(1,2))    # LSA
	y <- abs(coefficients(s)[-1,3])[1:(nProj/2)]                                
	plot( 
		x <- rep(1:length(y)), y, 
		xlab="Wine Regr, LSA variables", ylab="|t|", main="")
		abline(h=-qnorm(.025/length(y)), col="gray", lty=3)
		lines(lowess(x,y), col="red")
	half.normal.plot(y, height=5)
reset()

#    CCA of bigram variables
left  <-   1        : (nProj/2)
right <- (nProj/2+1):  nProj
ccw <- cancor(x.bigram.[,left], x.bigram.[,right]); 
plot(ccw$cor, xlab="Canonical Variable", ylab="Canonical Correlation")           # simccab.pdf

cl <- x.bigram.[, left] %*% ccw$xcoef		# canonical vars
cr <- x.bigram.[,right] %*% ccw$ycoef

regr.cl <- lm(rating ~ cl)
s.cl <- summary(regr.cl); s.cl

par(mfrow=c(1,2))    # Left bigram, after CCA
	y <- abs(coefficients(s.cl)[-1,3])[1:(nProj/2)]                                
	plot( 
		x <- rep(1:length(y)), y, 
		xlab="Wine Regr, Bigram variables (left, after CCA)", ylab="|t|", main="")
		abline(h=-qnorm(.025/length(y)), col="gray", lty=3)
		lines(lowess(x,y), col="red")
	half.normal.plot(y, height=5)
reset()

regr <- lm(rating ~ cl + cr)
s <- summary(regr); s









