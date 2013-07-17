### Eigenwords

UV <- read.csv("/Users/bob/Desktop/svd.csv"); dim(UV)

cols <- c("U7","U8")
plot(UV[,cols]); text(UV[,cols[1]],UV[,cols[2]], UV[,1])

### Regression

Regr <- read.table("/Users/bob/Desktop/regr.txt", header=TRUE); dim(Regr)
colnames(Regr)

# --- log first column
hist(log(Regr[,1]));
Regr[,1] <- log(Regr[,1])

# --- remove outlying sqft if present
hist(Regr[,"SqFt"]); which(Regr[,"SqFt"]>30000)

# --- scattplots
pairs(Regr[1:8])
pairs(Regr[c(1,9:15)])
pairs(Regr[c(1,16:22)])

# --- stepwise 
library(lars)

# --- standardized fit
step.fit  <- lars(as.matrix(Regr[,2:ncol(Regr)]), Regr[,1], type="stepwise", max.steps=10)
summary(step.fit)
cols <- as.numeric(step.fit$actions)
y <- Regr[,1]
X <- as.matrix(Regr[,cols])
fit <- lm(y~X)
summary(fit); plot(fit)


# --- other methods are similar to each other, but not to stepwise
lasso.fit <- lars(as.matrix(Regr[,2:ncol(Regr)]), Regr[,1], type="lasso"   , max.steps=15)

lars.fit  <- lars(as.matrix(Regr[,2:ncol(Regr)]), Regr[,1], type="lar"     , max.steps=15)

