
#base <- "/home/jason/genome_sequencing/StrainSeq_8.14.2019/"
base <- commandArgs(trailingOnly=TRUE) ##this accepts the dir input from terminal. trailing slash should be provided by terminal input
setwd(base)

#get fileList
dirs <- list.dirs(base, recursive=FALSE)
for (dir in dirs) {
  fList <- list.files(dir, pattern="*tsv", recursive = FALSE) #get fileList (should be 2 files)
  if (length(fList)==0) {
    dirs <- dirs[-which(dirs==dir)] #remove this dir
    next
  }
  tmp <- read.delim(paste0(dir,"/",fList[grep("bam_summary", fList)]), header=FALSE, stringsAsFactors = FALSE) #read in the main summary file
  if (dir==dirs[1]) { #set up dataframe if first iteration of loop
    data <- data.frame(Feature=tmp[,1], stringsAsFactors = FALSE)
  }
  
  #now get the rDNA info
  tmp2 <- read.delim(paste0(dir,"/",fList[grep("rDNA", fList)]), header=FALSE, stringsAsFactors = FALSE) #read in the rDNA file
  rVal <- as.numeric(tmp2[1,2]) #rDNA # of sequences
  
  ##now subtract this value from the Chr12 data
  tmp[which(tmp[,1]=="Chr12"),3] <- tmp[which(tmp[,1]=="Chr12"),3]-rVal
  
  data[,ncol(data)+1] <- tmp[,3]/tmp[,2] #calculate the ratio, append new column
  data[nrow(data),ncol(data)] <- 2*(rVal/17356)/data[which(data[,1]=="Chr12"),ncol(data)] #normalize to length then divide by the new average I calculated for chr12. Apply multiply by 2 correction
  
  
  
}

#inelegant way to extract the directory names
names <- strsplit(dirs,"/")
names2 <- vector()
for (i in 1:length(names)) {
  names2 <- c(names2,names[[i]][length(names[[i]])])
}

#rename things
colnames(data) <- c("Feature", names2)
data[nrow(data),1] <- "rDNA"

#normalize everything to Chr1. in hindsight this would all have been easier by working with a matrix type
#do not include the rDNA
data[-nrow(data),2:ncol(data)] <- apply(data[-nrow(data),2:ncol(data)], 2, function(x) {x/x[1]})

data[,2:ncol(data)] <- round(data[2:ncol(data)],digits=2)

write.table(data, file="copyNumberSummary.tsv", quote=FALSE, col.names=TRUE, row.names=FALSE, sep="\t")


