##################################################################################
# 
# Function definitions
#
##################################################################################

library(car)
library(xtable)    # latex tables
library(MASS)
library(leaps)

predictive.r2 <- function(regr) {  # returns all three
	n <- length(regr$residuals)
	aov <- anova(regr)
	ss <- aov$"Sum Sq"; 
	rss <- ss[length(ss)]; regr.ss <- sum(ss[-length(ss)]); tss <- regr.ss + rss
	k <- n - regr$df.residual  # k includes constant
	c(1-rss/tss,  1 - (rss/(n - k))/(tss/(n-1)), 1 - (rss/(n - 2*k))/(tss/(n-1)))
	}
	
	
rdirichlet <- function(a) {
    y <- rgamma(length(a), a, 1)
    return(y / sum(y))
	}

reset <- function() {
	# par(mfrow=c(1,1), mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)      # default
	par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,2.5,2,1)+0.1)  # bottom left top right
	}

coef.summary.plot <- function(sr, xlab, omit=1) { 				# omit intercept
	par(mfrow=c(1,2))  
		y <- abs(sr$coefficients[-(1:omit),3])
		x <- 1:length(y)     
		threshold <- -qnorm(.025/length(y))
		cat("Bonferroni at ",threshold,"\n");
		plot(x,y, xlab=xlab, ylab="|t|", main="", col="darkgray", cex=0.5)
		abline(h=threshold, col="black", lty=4)
		abline(h=sqrt(2/pi), col="black")
		lines(lowess(x,y,f=0.3), col="red")
		half.normal.plot(y)
	reset()
	}


half.normal.plot <- function(y, height=2) {
	k <- length(y)          
	plot(x <- qnorm(.5+(0:(k-1))/(2*k)), y <- sort(y), cex=0.5, col="darkgray",
		xlab="Normal Quantile", ylab="Sorted |t|"); 
	k <- floor(0.25*k)
	cat("Using lower 20% of cases (",k,")\n"); 
	regr <- lm(y[1:k] ~ x[1:k])
	abline(0,1, col="black", lty=2)
	abline(regr,col="red")
	text(0.25,height,paste("b =",round(coefficients(regr)[2],1)), cex=0.7)
	cat("t-stat and p-value are ", coefficients(summary(regr))[2,],"\n")
	}
	
jitter <- function(x) { x + 0.05 * sd(x) * rnorm(length(x)) }

correctly.ordered <- function(y, y.hat, n) { # percentage of random pairs correctly ordered.
	ii <- sample(1:length(y),n,replace=FALSE)
	jj <- sample(1:length(y),n,replace=FALSE)
	dy <- y[ii]-y[jj]
	df <- y.hat[ii]-y.hat[jj]
	plot(df,dy,xlab="Difference in Predicted Values", ylab="Difference in Actual Y");
	abline(a=0,b=1,col="red")
	abline(h=0, col="gray"); abline(v=0, col="gray")
	sum((dy*df)>0)
	}
	

################################################################
#
#   Cross validation variations
#
################################################################


cross.validate.mse <- function(data, n.folds=10,seed=23479,n.reps=1) { 	
	n <- length(data$y);
	i <- rep(1:n.folds,ceiling(n/n.folds))
	if (length(i) != n) cat("Note: n is not multiple of # folds.\n");
	i <- i[1:n]
	set.seed(seed)
	mse <- rep(0,n.folds*n.reps)
	for(kk in 1:n.reps) {
		ii <- sample(i, n)   # permute fold indices
		for(fold in 1:n.folds) {
			cat(fold,"\n");
			train <- (fold != ii)
			r <- lm(y ~ ., data=data[train,])
			test  <- (fold == ii)
			err <- data$y[test] - predict(r, newdata=data[test,]);
			mse[fold+(kk-1)*n.folds] <- sum(err^2)/sum(test)
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

