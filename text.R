##################################################################################
# Analysis of regressors 
##################################################################################

nProj <- 200

city  <- "Chicago"
file  <- paste("/Users/bob/C/text/text_src/temp/",city,"_bigram_regr.txt",sep="")

Data <- read.table(file, header=TRUE); dim(Data)

# --- check two code versions (use doboth in makefile; force same seeds prior to rand projection)
# X <- as.matrix(Data[,209:308])  # as computed within regressor.cc
# Y <- as.matrix(Data[,309:408])  #                    lsa.cc
# plot(cancor(X,Y)$cor)           # == 1

# --- about 87% have 100 or fewer tokens (many of which are punctuation)
nTokens  <- as.numeric(Data[,"n"])
fivenum(nTokens); quantile(nTokens,0.87)
hist(log10(nTokens))

# --- prices are very skewed, even after dropping many
price <- Data[,"Y"]
hist(price, breaks=50); 
hist(sort(price, decreasing=TRUE)[-(1:100)], breaks=50)
plot(price ~ nTokens)  # dominated by outliers

logPrice <- as.numeric(log(Data[,"Y"]))
hist(logPrice, breaks=30)
plot(logPrice ~ nTokens);  # some common lengths (vertical stripe)
   lines(lowess(nTokens,logPrice,f=0.1),col="red")
plot(logPrice ~ I(log(nTokens) )) 
   
sqft <- Data[,"SqFt"]  # too many missing; need log scale
plot(logPrice ~       sqft   )
plot(logPrice ~ I(log(sqft)) )  # clear coding error... min=1  "1sfam" in source)


parse.names<-c("n","SqFt","SqFt_Obs","Bedrooms","Bedroom_Obs","Bathrooms","Bathroom_Obs")

x.parsed. <- as.matrix(Data[,parse.names])
	x.parsed.[,"SqFt"] <- log(x.parsed.[,"SqFt"])   # transform to log (tiny improvement)
	colnames(x.parsed.)[2] <- "LogSqFt"

x.lsa.    <- as.matrix(Data[,paste("L",0:(nProj/2-1), sep="")])
x.bigram. <- as.matrix(Data[,paste("H",0:(nProj  -1), sep="")])

summary(regr.parsed        <- lm(logPrice ~ x.parsed.))
summary(regr.lsa           <- lm(logPrice ~ x.lsa.   ))
summary(regr.bigram        <- lm(logPrice ~ x.bigram.))

summary(regr.parsed.lsa    <- lm(logPrice ~ x.parsed. + x.lsa.))
summary(regr.parsed.bigram <- lm(logPrice ~ x.parsed. + x.bigram.))

summary(regr.all           <- lm(logPrice ~ x.parsed. + x.lsa. + x.bigram.))

# --- compare nested models
# adding lsa/bigram to parsed + bigram/lsa
anova(regr.parsed.lsa   , regr.all)
anova(regr.parsed.bigram, regr.all)

# adding lsa/bigram to parsed
anova(regr.parsed, regr.parsed.lsa)
anova(regr.parsed, regr.parsed.bigram)

# adding parsed to bigram/lsa
anova(regr.bigram, regr.parsed.bigram)
anova(regr.lsa   , regr.parsed.lsa)

# --- estimated parsed variables in place of originals
x.parsed.hat. <- x.parsed.
for(j in 1:ncol(x.parsed.)) {
	regr <- lm(x.parsed.[,j] ~ x.lsa.+ x.bigram.)
	cat("For variable", colnames(x.parsed.)[j], " R2 = ", summary(regr)$r.squared,"\n")
	x.parsed.hat.[,j] <- fitted.values(regr)
	}

#     just one at a time, skip n since redundant
colnames(x.parsed.)

#     n is worse when estimated, but all of the others improve
j <- 6;
colnames(x.parsed.)[j]
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
# Simulate data from topic model
##################################################################################

# --- functions

rdirichlet <- function(a) {
    y <- rgamma(length(a), a, 1)
    return(y / sum(y))
}

word.indices <- function(topics, P) {
	M <- ncol(P); m<-length(topics)
	z<-rep(0,m); 
	for(t in 1:m) z[t]<-sample(M,1,prob=P[topics[t],])
	z
}

# --- constants

beta <- c(1, 1, 1, 2, 2, 2, 3, 3,-2,-3) 	# topic weights
K    <- length(beta)		# number of topics

M <- 1000					# word types in vocabulary
n <- 4000					# num of documents


# --- topic distributions over words in K x M matrix P
alpha <- rep(.05, M)			# dirichlet parms, small alpha are spiky
P <- matrix(0, nrow=K, ncol=M)
for(k in 1:K) P[k,] <- rdirichlet(alpha)  

#     check dependence/correlation among distributions
round(P[1:10,1:12],3)
udv <- svd(P); plot(udv$d); round(cor(t(P)),2)
round(udv$u[,1],4); plot(udv$v[,1])


# --- X[i,] is distribution of topics in document i (random Dirichlet)
X <- matrix(0, nrow=n, ncol=K);
for(i in 1:n) X[i,] <- rdirichlet( rep(1/10,K) )

# Y is the response
Y <- X %*% beta + rnorm(n,sd=0.5)

# check "true" model; sum of the X's = 1 so no intercept, R2 around 90%
summary(regr <- lm(Y ~ X - 1))


# --- simulate documents; lengths of the documents (neg bin)
vocab <- paste("w",1:M,sep="")
m <- rpois(n,lambda=rgamma(n,30,1))

# words in doc assigned to topic, then pick word for that topic
docs <- rep("",n)
for(i in 1:n) {
	topics <- sample(K, m[i], replace=TRUE, prob=X[i,]); 
	docs[i] <- paste(round(Y[i],3), paste(vocab[word.indices(topics,P)], collapse=" "), sep=" ")
	}
	
# write words to file with response
write(docs,"/Users/bob/C/text/text_src/sim/sim.txt")


# --- Analysis of C++ results

# --- look at type frequencies, zipf plot
type.cts <- sort(scan("/Users/bob/Desktop/type_freq.txt"), decreasing=TRUE)

x <- log(1:length(type.cts)); y<-log(type.cts)
plot(x,y, xlab="log rank", ylab="log frequency")
abline(regr<-lm(y~x),col="red"); coefficients(regr)


# --- get regression data
file  <- "/Users/bob/C/text/text_src/temp/sim_regr.txt"

Data <- read.table(file, header=TRUE); dim(Data)

nTokens  <- as.numeric(Data[,"n"])		# matches m above
resp <- Data[,"Y"]						# matches Y above

nProj <- 50

x.lsa.    <- as.matrix(Data[,paste("L",0:(nProj/2-1), sep="")])
x.bigram. <- as.matrix(Data[,paste("H",0:(nProj  -1), sep="")])

summary(regr.lsa           <- lm(Y ~ x.lsa.   ))
summary(regr.bigram        <- lm(Y ~ x.bigram.))        # better fit
cor(fitted.values(regr.lsa),fitted.values(regr.bigram)) #   high corr

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
ccw <- cancor(x.bigram.[,left], x.bigram.[,right]); plot(ccw$cor)

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





