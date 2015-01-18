## build data from makefile  (make embedded_data.txt); make controls n lines

source("missing_data.R")
source("~/C/text/functions.R")

## read data file

rawData <- read.delim("embedded_data.txt", header=T, as.is=T); dim(rawData)

## distribution over prepositions, then subset for specific prepositions

(prep.counts <- table(rawData[,1]))

Subset <- rawData[ (rawData[,1] =="to") | (rawData[,1] == "of"),]; dim(Subset)

Subset[,1] <-  as.numeric(0 + (Subset[,1] == "to")); mean(Subset[,1])

## patch the columns with missing data ... just one indicator for each eigen set
## some columns have about 80% missing

Data <-  fill.missing.numerical(Subset); dim(Data)

missing.cols <- (1+ncol(Subset)):ncol(Data)
plot( m <- colMeans(Data[,missing.cols]), main="Missing Shares")
text(1:length(missing.cols), m, colnames(Data)[missing.cols])


## check conditioning; can see (log-log) some near singularity
## due to the missing data columns near the end (FFword,FP)
## (PNP,PNword)
udv <- svd(Data[,-c(1,82,85)]);
plot(udv$d, log="xy", main="Check for Singular")

V <-  udv$v; dim(V)

plot(V[,ncol(V)])


## use linear regression to screen out collinear terms in
## both train and test segments
train <- 1 == rbinom(nrow(Data),size=1,prob=0.5); test <- !train

summary(regr <- lm(Y ~ ., data=Data[train,]))
omit.train <- which(is.na(coef(regr)))

summary(regr <- lm(Y ~ ., data=Data[test,]))
omit.test <- which(is.na(coef(regr)))

omit <- union(omit.train,omit.test)

summary(ols.regr <- lm(Y ~ ., data=Data[train,-omit]))

calibration.plot(fitted(ols.regr), Data[train,"Y"])  # not well calibrated


## whole world regression; logistic does not converge!
##  train

summary(log.regr <- glm(Y ~ ., data=Data[train,-omit], family=binomial(logit)))

calibration.plot(fitted(log.regr), Data[train,"Y"])  # much better

tab <- table(as.factor(Data[train,"Y"]), as.factor(0.5 <fitted(log.regr)))

prop.table(tab)


##  test... logistic pred breaks

pred <- predict(log.regr,type="response",newdata=Data[test,-omit])  # blows up

# pred <- pmin(10,pmax(-10,pred))

calibration.plot(pred, Data[test,"Y"])

tab <- table(as.factor(Data[test,"Y"]), as.factor(0.5 < pred))
prop.table(tab)


##  ROC
library(pROC)

roc.train <- roc(Data[train,"Y"], fitted(log.regr), plot=T, col="black")
roc.test <- roc(Data[test,"Y"], pred, plot=T, add=T       , col="blue")

roc.ols.train <- roc(Data[train,"Y"], fitted(ols.regr),                           plot=T, add=T, col="red")
roc.ols.test  <- roc(Data[test ,"Y"], predict(ols.regr,newdata=Data[test,-omit]), plot=T, add=T, col="pink")

auc(roc.train)
auc(roc.test)
auc(roc.ols.train)
auc(roc.ols.test)


### -------------------------  automated  ----------------------------------
library(pROC)

split.half.indicator <- function(n) {
    half <- floor(n/2);
    b <- c(rep(T,half),rep(F,half+n-2*half))
    return(sample(b,n,replace=F))
}

auc.classifier <- function(prep0,prep1,data, seed=12736, plot.auc=F, balance=F) {
    set.seed(seed);
    data <- data[ (data[,"Y"] ==prep0) | (data[,1] == prep1),];
    data[,"Y"] <-  as.numeric(0 + (data[,"Y"] == prep1));
    data <- fill.missing.numerical(data);
    cat("Analysis data set has dimension ",dim(data),"\n");
    if(!balance) {
        test <- !(train <- 1 == rbinom(nrow(data),size=1,prob=0.5));
    } else {
        n0 <- length(i0 <- which(data[,"Y"]==0));
        n1 <- length(i1 <- which(data[,"Y"]==1));
        n <- min(n0,n1);
        cat("Balanced samples use ",n," items from each group.\n")
        i0 <- sample(i0,n,replace=F)
        i1 <- sample(i1,n,replace=F)
        b  <- split.half.indicator(n);
        train <- c(i0[b],i1[b]); test <- c(i0[!b],i1[!b])
    }
    omit.train <- which(is.na(coef(lm(Y ~ ., data=data[train,])))) # remove singular X
    omit.test  <- which(is.na(coef(lm(Y ~ ., data=data[test, ]))))
    omit <- union(omit.train,omit.test)
    cat("Fitting logistic regression omitting ", length(omit), "columns.\n")
    log.regr <- glm(Y ~ ., data=data[train,-omit], family=binomial(logit))
    pred <- predict(log.regr,type="response",newdata=data[test,-omit])
    roc.train <- roc(data[train,"Y"], fitted(log.regr), plot=plot.auc, col="black")
    roc.test <- roc(data[test,"Y"], pred, plot=F, add=plot.auc       , col="red")
    return(c(auc(roc.train), auc(roc.test), dim(data)))
}

auc.classifier("of","to",rawData, balance=T)  # test case, should be about 0.88

doit <- function(data, threshold=10000, seed=12672) {
    prep.counts <- table(data[,"Y"]);
    prep.counts <- prep.counts[threshold <= prep.counts]
    names <-  names(prep.counts)
    results <- diag(prep.counts)
    rownames(results) <- colnames(results) <- names
    for(i in 1:(length(names)-1)) {
        for(j in (i+1):length(names)) {
            x <- auc.classifier(names[i],names[j],data=data,seed=seed,balance=T)
            results[i,j]<-x[1]; results[j,i]<-x[2]
        }
    }
    return(results);
}

auc.table <- doit(rawData)




##  300000  16 Jan 2014, balanced
            at        by       for      from        in        of        on        to      with
at       17798     0.811     0.758     0.589     0.714     0.873     0.728     0.802     0.767
by       0.805     10754     0.604     0.769     0.765     0.910     0.783     0.709     0.791
for      0.758     0.603     23446     0.750     0.726     0.833     0.755     0.776     0.682
from     0.590     0.753     0.745     13310     0.699     0.841     0.696     0.672     0.586
in       0.707     0.746     0.727     0.681     51174     0.832     0.662     0.783     0.744
of       0.871     0.897     0.828     0.834     0.831     74069     0.856     0.878     0.850
on       0.726     0.779     0.747     0.683     0.663     0.857     21193     0.607     0.626
to       0.790     0.702     0.766     0.663     0.774     0.873     0.605     27086     0.765
with     0.764     0.773     0.680     0.579     0.741     0.851     0.626     0.758     21219


##  300000  15 Jan 2014
            at        by       for      from        in        of        on        to      with
at       17798     0.808     0.599     0.743     0.710     0.858     0.729     0.797     0.766
by       0.801     10754     0.784     0.777     0.762     0.770     0.782     0.657     0.785
for      0.597     0.773     23446     0.753     0.727     0.819     0.753     0.772     0.612
from     0.736     0.762     0.745     13310     0.697     0.830     0.541     0.744     0.754
in       0.706     0.756     0.724     0.684     51174     0.830     0.665     0.771     0.746
of       0.859     0.771     0.818     0.821     0.832     74069     0.694     0.866     0.836
on       0.722     0.776     0.750     0.542     0.660     0.695     21193     0.770     0.739
to       0.795     0.653     0.769     0.730     0.771     0.865     0.765     27086     0.765
with     0.763     0.778     0.609     0.750     0.738     0.831     0.737     0.760     21219
