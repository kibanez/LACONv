ReadBig <- function( File, sep="\t", header=TRUE, comment.char = "", skip = 0){

   tab5rows <- read.table( File, header = header, nrows = 5, stringsAsFactors=FALSE, skip=skip, sep=sep)
   classes <- sapply( tab5rows, class )
   Nrows <- system( paste("wc -l", File, sep=" "), intern=TRUE )
   Nrows <- as.numeric( unlist( strsplit( Nrows, split=" " ) )[ 1 ] )
   read.table( File, header = header,
   			  sep = sep,
              colClasses = classes,
              nrows = Nrows,
              comment.char = comment.char,
              skip = skip,
              , stringsAsFactors=FALSE)
}



# Parameters Parsing (from Command Line)
options <- commandArgs(trailingOnly = T)

cov.file <- as.character(options[1])
gc.file <- as.character(options[2])
gender <- as.character(options[3])

if (!file.exists(cov.file)){
 	write("control_coverage R function: the path to the covFile is not correct. A coverage outcome is required", stderr())
}

if (!file.exists(gc.file)){
 	write("control_coverage R function: the path to the gc.file is not correct. A GC content outcome is required", stderr())
}



tryCatch(
     {
       gc <- ReadBig(gc.file,header=FALSE)
     },

     error=function(e) {
         print(paste(e$message,gc.file,sep=" "))
         return(2)
     },

     warning=function(cond) {
           message(cond)
     },
     finally={
	  gc <- ReadBig(gc.file,header=FALSE)
    }
)


tryCatch(
     {
       cov <- ReadBig(cov.file,header=FALSE)
     },

     error=function(e) {
         print(paste(e$message,cov.file,sep=" "))
         return(2)
     },

     warning=function(cond) {
           message(cond)
     },
     finally={
          cov <- ReadBig(cov.file,header=FALSE)
    }
)
names(cov) <- cov[1,]
cov <- cov[-1,]

gc.content <- as.numeric(gc$V2)
gc.content <- round(gc.content,2)

# we compute each file independently

df <- data.frame(cov$Target,cov$total_coverage,cov$average_coverage)
df2 <- data.frame(cov$Target,cov$total_coverage,cov$average_coverage)

# average
num.reads.norm.avg <- c()
num.reads.singc.avg <- c()

for (i in 1:length(cov$Target)){ 
	index.gc <- which(gc.content == gc.content[i])
	# median of reads in windows with = GCcontent
	sum.avg <- c()
	sum.avg <- cov[,3][index.gc]
#	for (k in 1:length(index.gc)){
#		sum.avg <- c(sum.avg,cov[,3][index.gc[k]])
#	}

	# ACTHUNG !! as.numeric (!!!!)
	sum.avg <- as.numeric(sum.avg)
	median.gc.i.avg <- median(sum.avg)

	if ((median.gc.i.avg == 0) || is.na(median.gc.i.avg)){
		num.reads.norm.avg[i] <- 0
	}else{
		num.reads.norm.avg[i] <- as.numeric(cov[,3][i]) * (median(as.numeric(cov[,3]))) / median.gc.i.avg
	}
	#    num.reads.singc.totalcov[i] <- num.reads.totalcov[i] / median(num.reads.totalcov)
} # for i


# normalization of the num reads in order to compare across samples with different levels of seq.coverage
# read_norm = reads / median (reads in all the sample)

  num.reads.singc.avg <- as.numeric(cov[,3]) / median(as.numeric(cov[,3]))
  
  

##### CAMBIO NUEVO PARA PODER COMPARAR CONTROLES CON TESTS#########################################################
  # normalization II
  # Normalizar cada sample para tener mean = 0, std.deviation = 1

  num.reads.singc.avg2 <- scale(num.reads.singc.avg,center=T,scale=T)

############################################################################	

  median.num.reads.norm.avg <- median(as.numeric(num.reads.norm.avg))
  num.reads.norm.comparation.avg <- num.reads.norm.avg / median.num.reads.norm.avg


##### CAMBIO NUEVO PARA PODER COMPARAR CONTROLES CON TESTS#########################################################
  # normalization II
  # Normalizar cada sample para tener mean = 0, std.deviation = 1

  num.reads.norm.comparation.avg2 <- scale(num.reads.norm.comparation.avg,center=T,scale=T)

############################################################################	

  df <- cbind(df,num.reads.norm.comparation.avg,num.reads.singc.avg)
  df2 <- cbind(df2,num.reads.norm.comparation.avg2,num.reads.singc.avg2)


# Normalization regarding the # reads, coverage...platform makes the previous for loop.
# NOW is fundamental to normalize each one (each column, that represents each sample), so as to have mean 0 and std deviation 1 ==> in the same distribution of values

names(df) <- c("Target","total_coverage","average_coverage","gc_avg","singc_avg")
names(df2) <- c("Target","total_coverage","average_coverage","gc_avg","singc_avg")

#cov.out2 <- as.data.frame(as.matrix(df2),stringsAsFactors=F)

aux = paste("_normalized",gender,sep="_")
aux2 = paste("_normalized2",gender,sep="_")
izena <- paste(cov.file,aux,sep="")
izena2 <- paste(cov.file,aux2,sep="")
write.table(df,file=izena,quote=F,col.names=T,row.names=F,sep="\t")
write.table(df2,file=izena2,quote=F,col.names=T,row.names=F,sep="\t")