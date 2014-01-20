##################################################################################
# 
# Function definitions
#
##################################################################################

library(car)
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
	
	
	
# --- run the following to initialize
reset()  # sets up plot

