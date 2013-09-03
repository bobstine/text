##################################################################################
# Preliminaries 
##################################################################################

reset <- function() {
	par(mfrow=c(1,1), mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)  # bottom left top right
	}
	


##################################################################################
# Analysis of text regressors 
##################################################################################

# --- look at type frequencies, zipf plot   (zipf.pdf)
type.cts <- sort(scan("/Users/bob/Desktop/type_freq.txt"), decreasing=TRUE)

x<-1:length(type.cts); y<-type.cts
zipf.data <- data.frame(list(x=x,y=y,lx=log(x),ly=log(y)))

plot(y~x, xlab="rank", ylab="frequency", log="xy", data=zipf.data)
common.words <- c(".", ",", "and", "-", "in")
text(0.9*x[1:5],0.7*y[1:5],common.words,cex=c(1,1,0.5,1,0.5))

regr<-lm(ly~lx, data=zipf.data[1:500,]); coefficients(regr)
lx <- log(x<-c(1,5000)); y <- exp(predict(regr, data.frame(lx=lx)))
lines(x,y,col="red")


# --- analysis of regression models
nProj <- 200

city  <- "Chicago"
file  <- paste("/Users/bob/C/text/text_src/temp/",city,"_bigram_regr.txt",sep="")

Data <- read.table(file, header=TRUE); dim(Data)


# --- analysis of prices
y <- Data[,"Y"]
par(mfrow=c(1,2))                                             # prices.pdf
	hist(log10(y), breaks=30, main=" ",xlab="log10(Price)")
	qqnorm(log10(y), ylab="log10(Price)"); abline(a=mean(log10(y)),b=sd(log10(y)))
reset()


# --- frequencies of missing data
sum(Data[,"SqFt_Obs"])/nrow(Data)
sum(Data[,"Bedroom_Obs"])/nrow(Data)
sum(Data[,"Bathroom_Obs"])/nrow(Data)


# --- check two code versions (use doboth in makefile; force same seeds prior to rand projection)
# X <- as.matrix(Data[,209:308])  # as computed within regressor.cc
# Y <- as.matrix(Data[,309:408])  #                    lsa.cc
# plot(cancor(X,Y)$cor)           # == 1


# --- length of documents; about 87% have 100 or fewer tokens (many of which are punctuation)
nTokens  <- as.numeric(Data[,"n"])

mean(nTokens); fivenum(nTokens); quantile(nTokens,0.87)
boxplot(nTokens, horizontal=TRUE, xlab="Lengths of Descriptions")   # boxplot.pdf
hist(log10(nTokens))


# --- simple models for log of prices
logPrice <- as.numeric(log(Data[,"Y"]))
plot(logPrice ~ I(log(nTokens) )) 
lines(lowess(log(nTokens), logPrice, f=.3), col="red")


sqft  <- Data[,"SqFt"]  # too many missing; need log scale
baths <- Data[,"Bathrooms"]
beds  <- Data[,"Bedrooms"]

# --- plots of the parsed explanatory variables and response
par(mfrow=c(2,2), mar=c(4,4,1,1), mgp=c(2,1,0))
	plot(logPrice ~ nTokens)
	  text(500,6, paste("r=",round(cor(logPrice,nTokens),2)))
	plot(logPrice ~ baths,xlab="Number Bathrooms" ) 
	  text(7,6, paste("r=",round(cor(logPrice,baths),2)))
	plot(logPrice ~ beds ,xlab="Number Bedrooms"  )  
	  text(7,6, paste("r=",round(cor(logPrice,beds),2)))
	plot(logPrice ~ I(log(sqft)),xlab="log(Sq Ft)" )  # clear coding error... min=1  "1sfam" in source)
	  text(1,6, paste("r=",round(cor(logPrice,log(sqft)),2)))
	obs <- which(1==Data[,"SqFt_Obs"]); cor(logPrice[obs],log(sqft[obs]))
par(mfrow=c(1,1))


parse.names<-c("n","SqFt","SqFt_Obs","Bedrooms","Bedroom_Obs","Bathrooms","Bathroom_Obs")

x.parsed. <- as.matrix(Data[,parse.names])
	x.parsed.[,"SqFt"] <- log(x.parsed.[,"SqFt"])   # transform to log (tiny improvement)
	colnames(x.parsed.)[2] <- "LogSqFt"

summary(regr.parsed        <- lm(logPrice ~ x.parsed.))


# --- SVD variables, D
x.lsa.    <- as.matrix(Data[,paste("D",0:(nProj/2-1), sep="")])

summary(regr.lsa           <- lm(logPrice ~ x.lsa.   ))    # pvalue_a
	plot(coefficients(summary(regr.lsa))[,4], xlab="Singular Vector", ylab="P-value", main="")


# --- SVD variables, B
x.bigram. <- as.matrix(cbind(	Data[,paste("BL",0:(nProj/2-1), sep="")],
									Data[,paste("BR",0:(nProj/2-1), sep="")]  ))

summary(regr.bigram        <- lm(logPrice ~ x.bigram.))
	plot(c(0,1:100,1:100),                                 # pvalue_b
		coefficients(summary(regr.bigram))[,4], xlab="Singular Vector", ylab="P-value", main="")

summary(regr.bigram        <- lm(logPrice ~ x.bigram.[,1:(nProj/2)]))


# --- CCA of bigram left/right decomp

left  <-   1        : (nProj/2)
right <- (nProj/2+1):  nProj
ccw <- cancor(x.bigram.[,left], x.bigram.[,right])
plot(ccw$cor, xlab="Combination", ylab="Canonical Correlation")          # cca.pdf

cx <- x.bigram.[, left] %*% ccw$xcoef		# canonical vars
cy <- x.bigram.[,right] %*% ccw$ycoef

cor(cx[,1],cy[,1])

summary(regr.bigram        <- lm(logPrice ~ cx + cy ))


# --- nested regression models

summary(regr.parsed        <- lm(logPrice ~ x.parsed.         ))

summary(regr.parsed.lsa    <- lm(logPrice ~ x.parsed. + x.lsa.))

summary(regr.parsed.bigram <- lm(logPrice ~ x.parsed. + x.bigram.))

summary(regr.all           <- lm(logPrice ~ x.parsed. + x.lsa. + x.bigram.))

summary(regr.text          <- lm(logPrice ~             x.lsa. + x.bigram.))


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

M <- 2000					# word types in vocabulary
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
plot(x,y, xlab="log rank", ylab="log frequency")        # simzipf.pdf
abline(regr<-lm(y~x),col="red"); coefficients(regr)


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





