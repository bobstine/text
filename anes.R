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

  Cheney$name <- "Cheney" ; Cheney$ncol <-ncol(Cheney)-3
  Cheney$InterviewerTranscript <- tolower(as.character(Cheney$InterviewerTranscript))

  Pelosi$name <- "Pelosi" ; Pelosi$ncol <- ncol(Pelosi)-3
  Pelosi$InterviewerTranscript <- tolower(as.character(Pelosi$InterviewerTranscript))

  Roberts$name <- "Roberts" ; Roberts$ncol <- ncol(Roberts)-3
  Roberts$InterviewerTranscript <- tolower(as.character(Roberts$InterviewerTranscript))
	
dfs <- list(brown=Brown, cheney=Cheney, pelosi=Pelosi, roberts=Roberts)
	
	# --- separate into words
par(mfrow=c(2,2))
	for (d in 1:4)
	{	df <- dfs[[d]]
		df$words <- strsplit(df$InterviewerTranscript," ")
		df$nTokens <- as.vector( sapply(df$words,length) )
		df$words[[which.max(df$nTokens)]]
		hist(df$nTokens, xlab="Number of Word Tokens", main=df$name[1], xlim=c(0,60), breaks=20)
		text(40,500, paste("mean = ", round(mean(df$nTokens),1)))
		cat("Median number of word tokens,", df$name[1], median(df$nTokens)," mean ", mean(df$nTokens), "\n")
	}
reset()
		
	# --- count number of codes; related to text length?

par(mfrow=c(2,2))
	for (d in 1:4)
	{	df <- dfs[[d]]
		df$codes  <- df[,paste("Code",0:(df$ncol[1]-1),sep="")]
		df$nCodes <- apply(df$codes,1,function(x) df$ncol[1]-sum(is.na(x)))
		hist(df$nCodes, xlab="Number of Codes",  xlim=c(0,10), main=paste("Number of Codes,", df$name[1])); 
		text(6,600, paste("mean =",round(mean(df$nCodes),1)))   # 1.63
	# counts <- table(df$nCodes); p <- counts["2"]/counts["1"]
	# lines(counts["1"]*p^(0:ncol))    # drops faster than geometric
	}
reset()


	plot(jitter(df$nCodes), jitter(df$nTokens), main=name,
		xlab="Number of Assigned Codes",ylab="Number of Word Tokens")
	r <- lm(nTokens ~ poly(nCodes,2), data=df); summary(r)
	i <- order(df$nCodes)
	lines(df$nCodes[i], fitted.values(r)[i], col="red")

	df$words[which(0 == df$nCodes)]  # no text for 2


	# --- build vocabulary

	all.counts <- table(c(df$words, recursive=TRUE))
	cat ("Size of vocabulary for", name, length(all.counts))
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

