### Handle missing numerical data columns in R
library(stringr)

fill.missing.categorical <- function(data) {
	for(j in 1:ncol(data)) {
		if(is.factor(data[,j])) { # cat("Checking col ",j," for missing\n")
			i <- is.na(data[,j])
			if(any(i)) {	v <- as.vector(data[,j]);
				        	v[i] <- "Missing"; data[,j]<-as.factor(v)
				        }
		}}
	data
	}

## --- This function fills NA in numerical cols by avg of those present
##     and tacks on extra columns for indicators of missing data. Returns
##     data frame with results.
##
##     BEWARE: this version is designed to handle the NLP eigenword blocks

add.missing.indicators <- function(name) { # true if last two chars are _0
    return("_0" == str_sub(name,-2))
}

fill.missing.numerical <- function(data) {
    missing.names <- rep("0",ncol(data))    # empty space to fill
    missing.cols  <- matrix(0,nrow=nrow(data),ncol=ncol(data))
    col <- 1;
    names <-  colnames(data);
    for(j in 1:ncol(data)) {
        if(!is.factor(data[,j])) { # cat("Checking col ",j," for missing\n")
            i <- is.na(data[,j])
            if(any(i)) {
                m <- mean(data[,j],na.rm=TRUE);
                data[i,j] <- m;
                if(add.missing.indicators(names[j])) {
                    missing.names[col] <- paste("Miss.",names(data)[j],sep="")
                    missing.cols[,col]  <- 0+i;
                    col <- col+1
                }
            }}}
    missing.names <- missing.names[1:(col-1)]
    missing.cols  <- data.frame(missing.cols[,1:(col-1)]);
    names(missing.cols)<-missing.names
    data.frame(data,missing.cols)
}


#############################################
example <- function() {
	# make complete data
	example.df <- data.frame( cbind(x1=1,x2=1:10,x3=rnorm(10)),
                          lab=sample("UVW",10,replace=TRUE),
                          fac=sample("ABC",10,replace=TRUE))
	example.df

	# add some missing values
	example.df$x1[5]<-NA
	example.df$x2[c(1,4,6)]<-NA
	example.df$fac[9:10] <- NA
	example.df

	# fill in the data in various ways
	fill.missing.categorical(example.df)
	fill.missing.numerical(example.df)
	fill.missing(example.df)
	}
