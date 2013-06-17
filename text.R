

#############################################
Marginal distributions of types and POS
#############################################

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


#############################################
Comparison of cluster analysis of random projection matrix
#############################################

proj <- read.table("/Users/bob/Desktop/kmeans_data.txt"); dim(proj)

norm <- function(x) { sqrt(sum(x*x)) }

norm(proj[  1,1:50]);  norm(proj[  1,51:100])
norm(proj[ 10,1:50]);  norm(proj[ 10,51:100])
norm(proj[100,1:50]); norm(proj[100,51:100])

km <- kmeans(proj,200,iter.max=30)





