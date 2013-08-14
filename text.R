##################################################################################
# Analysis of regressors 
##################################################################################

nProj <- 200

city  <- "Chicago"
file  <- paste("/Users/bob/C/text/text_src/temp/",city,"_data.txt",sep="")

Data <- read.table(file, header=TRUE, row.names="Context"); dim(Data)


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

x.lsa.    <- as.matrix(Data[,paste("R",0:(nProj/2-1), sep="")])
x.bigram. <- as.matrix(Data[,paste("X",0:(nProj  -1), sep="")])

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

# --- canonical correlations are quite large, drop slowly
cc <- cancor(x.lsa., x.bigram.)
plot(cc$cor)

# --- cca of the left/right bigram variables; inverted hockey stick
#     Related to how many to keep?  
#     Clearly useful to trim left/right sets down:
#        only keep 'right' subspace that's not redundant with left
left  <-   1        : (nProj/2)
right <- (nProj/2+1):  nProj
cc <- cancor(x.bigram.[,left], x.bigram.[,right])
plot(cc$cor)



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





