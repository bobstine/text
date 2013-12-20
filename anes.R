#
#  Analysis of ANES data from LDC
#

dataPath <- "/Users/bob/C/text/text_src/anes/"
add.path <- function(file) { paste(dataPath,file,sep="") }

reset <- function() {
	# par(mfrow=c(1,1), mgp=c(3,1,0), mar=c(5,4,4,2)+0.1)      # default
	par(mfrow=c(1,1), mgp=c(1.5,0.5,0), mar=c(3,2.5,2,1)+0.1)  # bottom left top right
	}

reset()


#####################################################################################
#
#    Office recognition for Gordon Brown
#
#####################################################################################


Brown   <- read.csv(add.path("brown.csv")); dim(Brown)       # 2098 x 10, but last was empty
Cheney  <- read.csv(add.path("cheney.csv")); dim(Cheney)     # 
Pelosi  <- read.csv(add.path("pelosi.csv")); dim(Pelosi)     # 2098 x 12
Roberts <- read.csv(add.path("roberts.csv")); dim(Roberts)     # 2098 x 12

# --- text is in second column; convert to lower-case text
  
Brown$name <- "Brown" ; Brown$ncol <- ncol(Brown)-3
Brown$InterviewerTranscript <- tolower(as.character(Brown$InterviewerTranscript))
Brown$words <- strsplit(Brown$InterviewerTranscript," ")
Brown$nTokens <- as.vector( sapply(Brown$words,length) )

Cheney$name <- "Cheney" ; Cheney$ncol <-ncol(Cheney)-3
Cheney$InterviewerTranscript <- tolower(as.character(Cheney$InterviewerTranscript))
Cheney$words <- strsplit(Cheney$InterviewerTranscript," ")
Cheney$nTokens <- as.vector( sapply(Cheney$words,length) )

Pelosi$name <- "Pelosi" ; Pelosi$ncol <- ncol(Pelosi)-3
Pelosi$InterviewerTranscript <- tolower(as.character(Pelosi$InterviewerTranscript))
Pelosi$words <- strsplit(Pelosi$InterviewerTranscript," ")
Pelosi$nTokens <- as.vector( sapply(Pelosi$words,length) )

Roberts$name <- "Roberts" ; Roberts$ncol <- ncol(Roberts)-3
Roberts$InterviewerTranscript <- tolower(as.character(Roberts$InterviewerTranscript))
Roberts$words <- strsplit(Roberts$InterviewerTranscript," ")
Roberts$nTokens <- as.vector( sapply(Roberts$words,length) )



dfs <- list(brown=Brown, cheney=Cheney, pelosi=Pelosi, roberts=Roberts)

	
	# --- separate into words
par(mfrow=c(2,2))
	for (d in 1:4)
	{	df <- dfs[[d]]
		df$words[[which.max(df$nTokens)]]
		hist(df$nTokens, xlab="Number of Word Tokens", main=df$name[1], xlim=c(0,60), breaks=20)
		text(40,500, paste("mean = ", round(mean(df$nTokens),1)))
		cat("Median number of word tokens,", df$name[1], median(df$nTokens)," mean ", mean(df$nTokens), "\n")
	}
reset()
		
# --- count number of codes; related to text length?

Brown$codes  <-   Brown[,paste("Code", 0:(Brown$ncol[1]-1 ),sep="")]
Cheney$codes <-  Cheney[,paste("Code", 0:(Cheney$ncol[1]-1),sep="")]
Pelosi$codes <-  Pelosi[,paste("Code", 0:(Pelosi$ncol[1]-1),sep="")]
Roberts$codes<- Roberts[,paste("Code", 0:(Roberts$ncol[1]-1),sep="")]

Brown$nCodes   <- apply(  Brown$codes, 1, function(x)   Brown$ncol[1]-sum(is.na(x)))
Cheney$nCodes  <- apply( Cheney$codes, 1, function(x)  Cheney$ncol[1]-sum(is.na(x)))
Pelosi$nCodes  <- apply( Pelosi$codes, 1, function(x)  Pelosi$ncol[1]-sum(is.na(x)))
Roberts$nCodes <- apply(Roberts$codes, 1, function(x) Roberts$ncol[1]-sum(is.na(x)))

dfs <- list(brown=Brown, cheney=Cheney, pelosi=Pelosi, roberts=Roberts)

par(mfrow=c(2,2))
	for (d in 1:4)
	{	df <- dfs[[d]]
		hist(df$nCodes, xlab="Number of Codes",  xlim=c(0,10), main=paste("Number of Codes,", df$name[1])); 
		text(6,600, paste("mean =",round(mean(df$nCodes),1)))   # 1.63
	# counts <- table(df$nCodes); p <- counts["2"]/counts["1"]
	# lines(counts["1"]*p^(0:ncol))    # drops faster than geometric
	}
reset()


par(mfrow=c(2,2))
	for (d in 1:4)
	{	df <- dfs[[d]]
		plot(jitter(df$nCodes), jitter(df$nTokens), main=df$name[1], ylim=c(0,60),
			xlab="Number of Assigned Codes",ylab="Number of Word Tokens")
		r <- lm(nTokens ~ poly(nCodes,2), data=df); summary(r)
		i <- order(df$nCodes)
		lines(df$nCodes[i], fitted.values(r)[i], col="red")
	}
reset()


	# --- build vocabulary; count number unique

	for (d in 1:4)
	{	df <- dfs[[d]];
		all.counts <- table(c(df$words, recursive=TRUE))
		cat ("Size of vocabulary for", df$name[1], length(all.counts),"\n")
	}
	
	
	
	all.counts <- sort(all.counts, decreasing=TRUE)
	all.counts[1:20]
	all.types <- names(all.counts)


	# --- remove rare types (oov); 680 for brown

	oov.types <- all.types[(all.counts==1)|(all.counts == 2)]
	length(oov.types)

	types <- all.types[all.counts>2]; length(types)


	# --- create document/type matrix W

	equal <- function(a,b) { as.numeric(a == b) }

	#     test it
	i <- 2;
	df$InterviewerTranscript[[i]]
	w <- (outer(df$words[[i]], types, FUN="equal"))
	b <- apply(w,2,sum)
	sum(b)

	W <- matrix(0, nrow=nrow(df), ncol=length(types))
	for(i in 1:nrow(Brown)) W[i,] <- apply(outer(as.vector(df$words[[i]]),types,FUN="=="),2,sum)

	#    difference in these counts are oov tokens
	apply(W,1,sum)[1:10]
	df$nTokens[1:10]

	#    tack on oov
	W <- cbind(W,df$nTokens-apply(W,1,sum))
	colnames(W)<-c(types,"OOV")

	#    check that lengths match total counts 10453
	sum(W)
	sum(df$nTokens)
}

# --- LSA , after frequency and entropy normalization

entropy <- function(p) {    # 0 to 1
	p.norm <- p[p>0]/sum(p)
	-sum(log(p.norm)*p.norm)/log(length(p)) }
	
	
L <- W;
for(i in 1:nrow(L)) L[i,] <- L[i,]/sum(L[i,])
for(j in 1:ncol(L)) L[,j] <- (1-entropy(L[,j])) * L[,j]

# find missing for(i in 1:nrow(L)) if(any(is.na(L[i,]))) cat(i,"\n")

udv <- svd(L)

plot(udv$d[-length(udv$d)], log="y",
	main="Singular Values, LSA (Brown Office Recognition)", xlab="Index", ylab="Singular Value")


# --- Cluster documents using kmeans (0,10=other; 1,11=dk)

k <- 2
clusters <- kmeans(udv$u[,1:100], centers = k)
clusters$size
100 * clusters$betweenss/clusters$totss

table(clusters$cluster,Brown$Code0)

Brown$InterviewerTranscript[clusters$cluster==2]
Brown[clusters$cluster==2,c(2,3)]

