source("/Users/bob/C/text/functions.R")

path <- "/Users/bob/C/text/text_src/temp/movie_ratings/"

add.path <- function(s) paste(path,s,sep="")


##################################################################################
#
#  Raw data investigation
#
##################################################################################

ds  <- scan("/data/movies/scale_data/Dennis+Schwartz/id.Dennis+Schwartz")
jb  <- scan("/data/movies/scale_data/James+Berardinelli/id.James+Berardinelli")
sr  <- scan("/data/movies/scale_data/Scott+Renshaw/id.Scott+Renshaw")
srh <- scan("/data/movies/scale_data/Steve+Rhodes/id.Steve+Rhodes")

intersect(ds,jb); intersect(ds,sr); intersect(ds,srh)
                  intersect(jb,sr); intersect(jb,srh)
                  					intersect(sr,srh)

##################################################################################
#  type counts, zipf
##################################################################################

# --- look at type frequencies, zipf plot   (zipf.pdf)
#     amazingly linear, with slope 1
#     60,255 word types before running tokenize script, 38,886 after (at min count 3)
#     With higher count threshold 10, have half as many 19288

# VOCB: Full vocabulary has 82631 types from input of 4384649 tokens on 5006 lines.                                     
# VOCB: Thresholded vocabulary of 38815 types with token count 4384649 tokens.                                                                                  
# VOCB: Position of OOV type is 9 with frequency 54342                                                                                                          
# MAIN: Vocabulary has 38815 types from 4384649 tokens, with 54342 OOV.  Most frequent are:                                                                      
#       ","->233729 "the"->213850 "."->184522 "a"->108441 "of"->97058 
#       "and"->94661 "to"->92188 "is"->72036 "in"->60992 "OOV"->54342                                      

# MAIN: regressor --file=regr_data_2.txt -min_frequency=10 --adjustment=b                                                                                                                                   
# VOCB: Full vocabulary has 55450 types from input of 2148164 tokens on 5006 lines.                                        
# VOCB: Thresholded vocabulary of 10931 types with token count 2148164 tokens.                                                                                            
# VOCB: Position of OOV type is 2 with frequency 104043                                                                                                                   
# MAIN: Vocabulary has 10931 types from 2148164 tokens, with 104043 OOV.  
# Most frequent 10 types are:                                                                     
#  "the"->112068 ","->108244 "OOV"->104043 "."->102043 "of"->55398 
#   "a"->55045 "and"->47942 "to"->43916 "is"->39547 "in"->28447                                      


type.cts <- sort(scan(add.path("type_freq.txt")), decreasing=TRUE)		

	x<-1:length(type.cts); y<-type.cts
	zipf.data <- data.frame(list(x=x,y=y,lx=log(x),ly=log(y)))

	plot(y~x, xlab="rank", ylab="frequency", log="xy", data=zipf.data)
	common.words <- c(",", "the",".","a","of","and", "to","is","in","OOV")
	text(0.9*x[1:5],0.7*y[1:5],common.words[1:5],cex=0.5)

	regr<-lm(ly~lx, data=zipf.data[1:500,]); coefficients(regr)
	lx <- log(x<-c(1,5000)); y <- exp(predict(regr, data.frame(lx=lx)))
	lines(x,y,col="red")


##################################################################################
#
# Response and document lengths
#
##################################################################################

# --- read data generated by 'regressor' (make dore)

	Data    <- read.table(add.path("lsa_ym.txt"), header=TRUE); dim(Data)

	n <- nrow(Data)
	logRating <- Data[,"Y"]    # file holds log prices
	rating    <- exp(Data[,"Y"])
	nTokens   <- Data[,"m"]
	logTokens <- log(nTokens)
	

# --- reviewer information
#     dropping 4 that are zero rated for Scott Renshaw  11961, 1391, 2790, 3285
	reviewer <- as.factor(c(rep("DS",1028),rep("JB",1308),rep("SRe",903-4),rep("SRh",1771)))
	col.r    <-  c(rep("red",1028),rep("seagreen",1308),rep("blue",903-4),rep("violet",1771))

# --- lengths (m) average was 876 from merged, down to 429 using subj files
	mean(nTokens); fivenum(nTokens)
	hist(log10(nTokens), breaks=20)

# --- analysis of ratings, little association with length  (r = -0.05)
#     no reason for log scale on ratings, and indeed skews the ratings
	hist(rating, breaks=20)
	hist(logRating, breaks=20);		
	plot(rating ~ nTokens,col=col.r ); cor(rating,nTokens)

	plot(rating ~ logTokens,col=col.r ); cor(rating,logTokens)  # a very large outlier
	rating.s <- rating[-450]; logTokens.s <- logTokens[-450]    # cor = 0.19
	plot(rating.s ~ logTokens.s,col=col.r ); cor(rating.s,logTokens.s)  # without
	
	
# --- reviewer explains about 4.4% of ratings, 3% of log ratings (not much reason to log)
#     write more about better movies, on average
#	  some interaction, but tiny improvement to fit
	boxplot(logRating~reviewer)
	summary(regr.a <- lm(rating~reviewer ))  # 4.3%
	summary(regr.b <- lm(rating~reviewer + logTokens))  # 11%  but only 4.7% with subj data

##################################################################################
#
#     LSA
#
##################################################################################

	file    <- add.path("lsa_cca_1000_p4.txt")
	LSA     <- as.matrix(read.table(file, header=TRUE)); dim(LSA)	# 5006 rows

# --- spectrum from random matrix is much less power-law with tf-idf filtering
	sv <- as.vector(scan(add.path("lsa_cca_1000_p4_d.txt")))
	plot(sv, log="xy", xlab="Component", ylab="Singular Value")

	# sv.all 	0.999999 0.295073 0.275427 0.238178 0.226597 0.216978
	# sv sort	0.764017 0.259918 0.233001 0.185635 0.179704 0.171391
	# sv subj   1.000000 0.344491 0.241269 0.233073 0.217952 0.215022
	
	plot(LSA[,1],col=col.r, cex=0.5)

	file    <- add.path("lsa_cca_1000_p4_v.txt")
	V       <- as.matrix(read.table(file, header=TRUE)); dim(V)	# 10,930 rows, 25 cols
	
	plot(-V[,1], log="xy", xlab="Word", cex=0.4)
	
# --- LSA analysis 		raw			log
#		wo reviewer		0.33		0.26
#       w  reviewer		0.344		0.275

	p      <- 500
	lsa    <- as.matrix(LSA[,1:p])
	
	sr.d <- summary(regr.d <- lm(rating ~ lsa)); sr.d
	predictive.r2(regr.d)  # 32% @ 100, 33% @ 200  but 47% with subj data!!!


	sr.e <- summary(regr.e <- lm(rating ~ reviewer + lsa)); sr.e
	predictive.r2(regr.e)  # 32% @ 100, 33% @ 200 and  48% with subj data
	
	# quartz(width=6.5,height=3); reset()
	coef.summary.plot(sr.d, "LSA Component", omit=6) 
	

# --- sequence of R2 statistics from C++  (watch for """ in C output)
	lsa.fit<- read.table(add.path("lsa_regr_fit_no_m_for.txt"),header=T)
	plot(lsa.fit[,"AICc"], cex=0.4, xlab="Model Size", ylab="AICc")
	plot(lsa.fit[,"adjR2"], cex=0.4, xlab="Model Size", ylab="Adjusted R2")
	
	# change names (legacy C++ labels (which are words) with types)
	# rownames(lsa.fit) <- c("tokens",paste("lsa",1:(nrow(lsa.fit)-1), sep=""))  
	rownames(lsa.fit) <- c(paste("lsa",1:(nrow(lsa.fit)), sep=""))  
	
	quartz(height=3.5, width=6); reset()
	plot(word.fit[,"AICc"], type="l", xlab="Features", ylab="AICc",    # [ aic.pdf portion ]
			lty=3, ylim=range(lsa.fit[,"AICc"]))
	lines(c(opt.k,opt.k), c(0,4300), col="gray")
	lines(lsa.fit[,"AICc"]) 
	lsa.fit[opt.lsa <- which.min(lsa.fit[,"AICc"]),]; opt.lsa
	lines(c(opt.lsa,opt.lsa), c(0,3250), col="gray")
	
	p <- 523;
	lsa    <- as.matrix(LSA[,1:p])
	sr.lsa <- summary(regr.lsa <- lm(logPrice ~ poly(nTokens,5) + lsa , x=TRUE, y=TRUE)); sr.lsa
	predictive.r2(regr.lsa)


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


