setwd("~/Рабочий стол/R data/resting-state dataset") 
library(data.table)

result_O <- read.table(file = "result_O.txt", header = T, sep=";")
#result_O <- data.frame(result_O[,-1], row.names=result_O[,1])
result_C <- read.table(file = "result_C.txt", header = T, sep = ";")


#get Raven data
raven <- read.table('all_twins.csv', header=TRUE, dec='.', sep=';')
names(raven)[names(raven)=="twinid"] <- "ID"

# regress mem out
raven <- na.omit(raven)
temp <- lm(pcrvtot~pccbtot, raven)
IQ_MEM <- residuals(temp)
raven$pure_IQ<-IQ_MEM
raven$pccbtot <- NULL

# моржируем
data_O <- merge(result_O, raven, by="ID")
data_C <- merge(result_C, raven, by="ID")

data_C <- data_C[!duplicated(data_C$ID),]
data_O <- data_O[!duplicated(data_O$ID),]

# split the twin data

data_O[["X"]] <- trunc(runif(nrow(data_O), 0, 2))
split <- split(data_O, data_O$X %% 2)
data_O <- split[[1]]

data_C[["X"]] <- trunc(runif(nrow(data_C), 0, 2))
split <- split(data_C, data_C$X %% 2)
data_C <- split[[1]]

data_O$X <- NULL
data_C$X <- NULL

# add GFP information
#GFP_alpha <- read.delim(file = "GFP_delta.txt", header = T, sep = " ", dec = " ")
#GFP_alpha <- read.table("GFP_alpha.txt", sep = ";", dec = '.', na.strings = c('???', 'NA'), header = TRUE)
#test <- readLines("GFP_alpha.txt")
#test1 <- strsplit(test, split = "\t")
#test2 <- data.frame(matrix(unlist(test1), nrow=191, byrow=T),stringsAsFactors=FALSE)
#test <- as.data.frame(do.call(rbind, strsplit(test, split = "")))

# сохраняем то, что получилось в папку реопзитория git
setwd("~/Рабочий стол/R data/Connectivity Project/Connectivity") 
data_O <- write.csv2(data_O, 'data_O.csv', row.names=F)
data_C <- write.csv2(data_C, 'data_C.csv', row.names = F)
