###arguments are: path/to/vcf.vcf, Q-score cutoff, setA (optional), setB (optional). The sets are numbers as appear in the vcf (1,2,3 vs. 4,5,6 for example)

#user-defined features
args2 <- commandArgs(trailingOnly=TRUE)

args <- c("",
          0,
          "",
          "",
          "/home/jason/annotations/SacCer3/SacCer3.fasta",
          "/home/jason/scripts/annotate-vcf/simplified-gene-coordinates_with_introns.txt")
    
if (length(args2)>0) {
    for (i in 1:length(args2)) {
      args[i] <- args2[i]
    }
  }

library(Biostrings)
getGeneticCode("3") -> GENETIC_CODE_M #yeast mito genetic code
wDir <- gregexpr(pattern="/", args[1])[[1]]
setwd(substr(args[1], start=1, stop=wDir[length(wDir)]))
fName <- args[1]
cutoff <- as.numeric(args[2]) #this is the quality score, which is universal among vcf files. I have some legacy lines of code in here that use the individual sample p-value, but only freebayes seems to output that

#now load in the vcf file. First load it in as a 1-column file to find the starting position of the matrix, then load it in again as a truncated file so it fits matrix format
vcfTMP <- read.delim(file=fName, header=TRUE, stringsAsFactors = FALSE)
vcf <- read.delim(file=fName, header=TRUE, stringsAsFactors = FALSE, skip=which(vcfTMP=="#CHROM"))
rm(vcfTMP)

gtf <- read.delim(file=args[6], header=TRUE, stringsAsFactors = FALSE)
fasta <- readDNAStringSet(filepath=args[5], format="fasta")

#now define columns
a <- c(10:(ncol(vcf))) #this vector specifies the columns containing genotypes. Re-use this later on.

#now define sets relative to a
setA <- a[as.numeric(strsplit(args[3],split=",")[[1]])]
setB <- a[as.numeric(strsplit(args[4],split=",")[[1]])]

getDifferences <- function(vcf, list1,list2) { #this function takes any two sets (each set can itself be a set) in terms of which columns, then finds everything that each set has in common, then finds things that are uncommon between the sets
  
  #first create a table identifying alleles that each list have in common
  if (length(list1) > 1) {
    set1 <- apply(vcf[,list1],1,function(x) {abs(max(x) - min(x)) == 0})
  } else {
    set1 <- rep(TRUE,nrow(vcf))
  }
  
  if (length(list2) > 1) {
    set2 <- apply(vcf[,list2],1,function(x) {abs(max(x) - min(x)) == 0})
  } else {
    set2 <- rep(TRUE,nrow(vcf))
  }
  
  #now, for each position that the lists share in common, retrieve the allele value
  #note that this retrieves the value from the first column of the list. This should be fine since I already subsetted on things they have in common
  alleles1 <- rep(NA,nrow(vcf))
  alleles1[which(set1==TRUE)] <- vcf[which(set1==TRUE),list1[1]]
  alleles2 <- rep(NA,nrow(vcf))
  alleles2[which(set2==TRUE)] <- vcf[which(set2==TRUE),list2[1]]
  
  return(which(alleles1!=alleles2))
}

getSimilarities <- function(vcf, list1,list2) { #this function takes any two sets (each set can itself be a set) in terms of which columns, then finds everything that each set has in common, then finds things that are also common between the sets
  
  #first create a table identifying alleles that each list have in common
  if (length(list1) > 1) {
    set1 <- apply(vcf[,list1],1,function(x) {abs(max(x) - min(x)) == 0})
  } else {
    set1 <- rep(TRUE,nrow(vcf))
  }
  
  if (length(list2) > 1) {
    set2 <- apply(vcf[,list2],1,function(x) {abs(max(x) - min(x)) == 0})
  } else {
    set2 <- rep(TRUE,nrow(vcf))
  }
  
  #now, for each position that the lists share in common, retrieve the allele value
  #note that this retrieves the value from the first column of the list. This should be fine since I already subsetted on things they have in common
  alleles1 <- rep(NA,nrow(vcf))
  alleles1[which(set1==TRUE)] <- vcf[which(set1==TRUE),list1[1]]
  alleles2 <- rep(NA,nrow(vcf))
  alleles2[which(set2==TRUE)] <- vcf[which(set2==TRUE),list2[1]]
  
  return(which(alleles1==alleles2 & alleles1!=0)) #note that these two statements are sufficient to ensure that neither set is simply the reference allele
}


findAssociatedGene <- function(position,chr, refVal, altVal) {
  #limitations: Does not deal with when a SNP matches to multiple genes (overlapping). Does not deal with multiple SNPs called for the same position in a single row (alt allele name will look like snpA,snpB)
  answers <- vector()
  #first subset data by chromosome. To speed things up could do this once ahead of itme
  #generate a 2-column table of SNP position relative to start and stop of all genes on chr
  position-gtf[gtf$Chr==chr,c(4,5)] -> tmp
  
  #find row which has the SNP lying within the gene intervals
  which(tmp[,1]>0 & tmp[,2]<0) -> row
  if (length(row)>0) {
    gtfRow <- as.numeric(row.names(tmp)[row]) #find row of matching gene
    answers[1] <- gtf[gtfRow,3] #report gene match
    #ensure that only a single SNP is listed as an alternative. If multiple are listed, break the script. I could alternatively just use the first value, but that might give misleading output
    if (length(grep(altVal, pattern=","))>0) {
      answers[2:5] <- "MULTIPLE"
      return(answers)
    }
    #also ensure that <DEL> or <INS> or other weird flags are not present
    if (length(grep(altVal, pattern=">"))>0) {
      answers[2:5] <- "LARGE_CHANGE"
      return(answers)
    }
    
    #extract whole sequence of the matched gene
    sequence <- fasta[[which(names(fasta)==chr)]][gtf[gtfRow,4]:gtf[gtfRow,5]]
    
    #get relative position of reference value
    relPos <- tmp[row,1]+1 #gets the relative position in the gene, 1 is added such that ATG starts at position 1, not 0
    relPos <- relPos:(relPos+(nchar(refVal)-1)) #gets the end of the reference position in case it's an indel
    
    #ensure that lookup value in fasta matches the reference value in the vcf
    fastaVal <- tryCatch({as.character(sequence[relPos])},
                         error=function(cond) {return("")})
    
    if (fastaVal==refVal) {
      #extract intron-corrected sequence coordinates
      coords <- as.numeric(strsplit(gtf[gtfRow,7],split="_")[[1]])
      coords2 <- c(coords[1]:coords[2])
      if (length(coords)>2) {
        for (coorbreaks in 2:(length(coords)/2)) {
          coords2 <- c(coords2, coords[(2*coorbreaks)-1]:coords[2*coorbreaks])
        }
      }
      #now get the intron-corrected sequence
      sequence2 <- fasta[[which(names(fasta)==chr)]][coords2]
      answers[2] <- length(sequence2)/3 #report protein length
      
      #get intron coordinates
      fullCoords <- gtf[gtfRow,4]:gtf[gtfRow,5]
      intronCoords <- fullCoords[-which(fullCoords %in% coords2)]
      
      #check whether the reference allele is fully contained within the intron. Use the relpos vector from earlier in cases of long allele variants
      relintronCoords <- intronCoords-min(fullCoords)+1
      if (any(relPos %in% relintronCoords)) {
        #check whether the splice junction is affected
        if (any(relPos %in% relintronCoords[1]) | any(relPos %in% relintronCoords[length(relintronCoords)])) {
          answers[3:5] <-"Splice Junction"
        } else {
          answers[3:5] <- "Intron"
        }
      } else {
        #create an alternative sequence. If the length changes, then I need to do it explicitly as below
        as.character(sequence) -> tmpSeq
        sequence3 <- DNAString(paste0(substr(tmpSeq,1,relPos[1]-1),altVal, substr(tmpSeq,start=relPos[length(relPos)]+1,stop=nchar(tmpSeq))))
        #now, at this point, sequence1 is the reference unspliced transcript, sequence2 is the reference spliced transcript, sequence3 is the alternative unspliced transcript
        #create sequence 4 which is the alternative transcript after splicing. This loop should not run if the alternative allele affects the intron at all
        if (length(relintronCoords)>0) {
          sequence4 <- sequence3[-relintronCoords]
        } else {
          sequence4 <- sequence3
        }
        #create a new relPos which has been adjusted for the intron removal
        #must be careful to adjust only for introns 5' of relPos, allowing gaps
        if (any(relPos>relintronCoords)) {
          relPos <- relPos-length(relintronCoords[which(relintronCoords<relPos)])
        }
        
        #declare genetic code for later
        if (chr=="mitochondrion") {
          code <- GENETIC_CODE_M
        } else {
          code <- GENETIC_CODE
        }
        
        #check for ATG to figure out complementarity
        if (countPattern("ATG", sequence2[1:3])==1) {
          answers[3] <- paste(as.character(sequence2[relPos]), min(relPos), altVal,sep="_")
          protA <- translate(sequence2, genetic.code=code)
          protB <- translate(sequence4, genetic.code=code)
          codons <- sort(unique(floor(((relPos)-1)/3)+1))
          codonsB <- codons
          if (max(codons)>length(protB)) {
            codonsB <- codons[codons<=length(protB)]
          }
          answers[4] <- paste(as.character(protA[codons]), min(codons), as.character(protB[codonsB]),sep="_")
          answers[5] <- if ((abs(nchar(altVal)-nchar(refVal)) %% 3)>0) {"Frameshift"
            } else if (as.character(protA[codons])==as.character(protB[codonsB])) {"Silent"
            } else {"Missense"
                }
        } else { 
          #RC everything and the rest should work
          sequence2 <- reverseComplement(sequence2)
          sequence4 <- reverseComplement(sequence4)
          relPos <- (nchar(sequence2)-relPos)+1
            
          answers[3] <- paste(as.character(sequence2[relPos]), min(relPos), as.character(reverseComplement(DNAString(altVal))),sep="_")
          protA <- translate(sequence2, genetic.code=code)
          protB <- translate(sequence4, genetic.code=code)
          codons <- sort(unique(floor(((relPos)-1)/3)+1))
          codonsB <- codons
          if (max(codons)>length(protB)) {
            codonsB <- codons[codons<=length(protB)]
          }
          answers[4] <- paste(as.character(protA[codons]), min(codons), as.character(protB[codonsB]),sep="_")
          answers[5] <- if (as.character(protA[codons])==as.character(protB[codonsB])) {"Silent"
            } else if ((abs(nchar(altVal)-nchar(refVal)) %% 3)>0) {"Frameshift"
            } else {"Missense"
            }
        }
      }
    } else {
      answers[2:5] <- "MISMATCH TO FASTA"
    }
  } else {
    tmp1 <- which(tmp[,2]==min(tmp[which(tmp[,2]>0),2])) #this returns the gene that the SNP is 3' of
    tmp2 <- which(tmp[,1]==max(tmp[which(tmp[,1]<0),1])) #this returns the upstream gene
    #now find the nearest downstream gene
    answers[1] <- paste(gtf[as.numeric(row.names(tmp)[tmp1]),3],tmp[tmp1,2], abs(tmp[tmp2,1]),gtf[as.numeric(row.names(tmp)[tmp2]),3],sep="_")
    answers[2:5] <- NA
  }
  return(answers)
}

#now add annotations and resave starting file
genes <- vector()
protlengths <- vector()
dnas <- vector()
prots <- vector()
types <- vector()

print("Processing SNPs...")
for (i in 1:nrow(vcf)) {
  if (i%%100==0) {
    print(paste0(i, " of ", nrow(vcf), " complete"))
  }
   answers <- findAssociatedGene(vcf[i,2],vcf[i,1], vcf[i,4], vcf[i,5])
   genes[i] <- answers[1]
   protlengths[i] <- answers[2]
   dnas[i] <- answers[3]
   prots[i] <- answers[4]
   types[i] <- answers[5]
}

vcfAnnotated <- cbind(vcf, GENE=genes, ProtLengths=protlengths, DNAMut=dnas, ProtMut=prots, Type=types)
name <- "Annotated"


#now reduce the vcf data to just a binary genotype. ASSUMES THE FIRST DELIMITER IN COLUMN 10 IS THE GENOTYPE. This has been true in all vcf files I've seen but it's still hardcoded
#vcf files which have a . in place of snps they didn't assess will get an NA value below. I change these to 0 for convenience but it's an assumption and is not as rigorous as assessing the SNPs directy in the SNP caller
vcf2 <- vcfAnnotated[1:(nrow(vcfAnnotated)-1),] #vcf2 will contain just the binary genotype calls and will be used for finding indices later
vcf2[,a] <- apply(vcf2[,a,drop=FALSE],c(1,2),function(x) {strsplit(x,":")[[1]][[1]]})
#now, convert the genotype strings to numbers
#haploids just need to be numeric, diploids need to be recoded as 0 (0/0), 1 (1/0 and 0/1), and 2 (1/1)
if (length(grep(pattern="/", x=vcf2[1,a]))>0) {
  vcf2[,a][vcf2[,a]=="0/0"] <- 0
  vcf2[,a][vcf2[,a]=="1/0"] <- 1
  vcf2[,a][vcf2[,a]=="0/1"] <- 1
  vcf2[,a][vcf2[,a]=="1/1"] <- 2
}
vcf2[,a] <- apply(vcf2[,a,drop=FALSE],c(1,2),function(x) {as.numeric(x)})
vcf2[,a][is.na(vcf2[,a])] <- 0


vcfAnnotated[1+nrow(vcfAnnotated),1] <- "Total SNPs"
vcfAnnotated[nrow(vcfAnnotated),a] <- apply(vcf2[,a], 2, function(x){length(which(x>0))})
write.table(vcfAnnotated, file=paste0(fName,"_",name,".txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(vcfAnnotated[which(as.numeric(vcfAnnotated[,6])>cutoff),], file=paste0(fName,"_",name, "_Q",cutoff, ".txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(vcfAnnotated[which(as.numeric(vcfAnnotated[,6])>cutoff & vcfAnnotated$Type!="Silent"& vcfAnnotated$Type!="Intron"),], file=paste0(fName,"_",name, "_Q",cutoff, "_Deleterious", ".txt"), sep="\t", quote = FALSE, row.names = FALSE)

write.table(vcf2, file=paste0(fName,"_",name,"_simple.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(vcf2[which(as.numeric(vcf2[,6])>cutoff),], file=paste0(fName,"_",name, "_Q",cutoff, "_simple.txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(vcf2[which(as.numeric(vcf2[,6])>cutoff & vcf2$Type!="Silent"& vcf2$Type!="Intron"),], file=paste0(fName,"_",name, "_Q",cutoff, "_Deleterious", "_simple.txt"), sep="\t", quote = FALSE, row.names = FALSE)



#now, create a summary of any SNPs that appear in only some of the strains tested, without respect to user-defined sets
name <- "All_shared"

#get indices of things that do not differ at all
indices <- which(apply(vcf2[,a],1,function(x) {abs(max(x) - min(x)) == 0})==TRUE)
unsharedVcf <- vcfAnnotated[indices,]
write.table(unsharedVcf, file=paste0(fName,"_",name,".txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(unsharedVcf[which(unsharedVcf[,6]>cutoff),], file=paste0(fName,"_",name, "_Q",cutoff, ".txt"), sep="\t", quote = FALSE, row.names = FALSE)



#now, create a summary of any SNPs that appear in all of the strains tested
name <- "All_not_shared"

#get indices of things that differ at all
indices <- which(apply(vcf2[,a],1,function(x) {abs(max(x) - min(x)) != 0})==TRUE)
unsharedVcf <- vcfAnnotated[indices,]
write.table(unsharedVcf, file=paste0(fName,"_",name,".txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(unsharedVcf[which(unsharedVcf[,6]>cutoff),], file=paste0(fName,"_",name, "_Q",cutoff, ".txt"), sep="\t", quote = FALSE, row.names = FALSE)



#legacy line2 of code in case I ever want to resurrect p-value cutoffs per genotype
#apply(unsharedVcf[,a[1],drop=FALSE],1,function(x) {any(as.numeric(strsplit(strsplit(x,":")[[1]][[8]],",")[[1]])<cutoff)}) -> truthsA
#now further subset the data
#unsharedVcf <- unsharedVcf[which(truthsA),]


#now, provide two sets and find any alleles that differ between the samples. This is relative to the sample order, not actual columns, to make it user-friendly.
name <- paste0("Set_",args[3],"_vs_",args[4],"_Differences")
getDifferences(vcf2,setA,setB) -> indices
subVcf <- vcfAnnotated[indices,]
write.table(subVcf, file=paste0(fName,"_",name,".txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(subVcf[which(subVcf[,6]>cutoff),], file=paste0(fName,"_",name, "_Q",cutoff, ".txt"), sep="\t", quote = FALSE, row.names = FALSE)

#use a massive one-liner to split the sample information for one column and 
#apply(subVcf[,setA[1],drop=FALSE],1,function(x) {any(as.numeric(strsplit(strsplit(x,":")[[1]][[8]],",")[[1]])<cutoff)}) -> truthsA
#apply(subVcf[,setB[1],drop=FALSE],1,function(x) {any(as.numeric(strsplit(strsplit(x,":")[[1]][[8]],",")[[1]])<cutoff)}) -> truthsB
#now further subset the data
#subVcf2 <- subVcf[which(truthsA & truthsB),]


#####now do it again but for things they have in common AND NOT REFERENCE ALLELES
#now, provide two sets and find any alleles that differ between the samples
#re-use setA and setB from above, but that could be modified here

name <- paste0("Set_",args[3],"_vs_",args[4],"_Shared")
getSimilarities(vcf2,setA,setB) -> indices
subVcf <- vcfAnnotated[indices,]
write.table(subVcf, file=paste0(fName,"_",name,".txt"), sep="\t", quote = FALSE, row.names = FALSE)
write.table(subVcf[which(subVcf[,6]>cutoff),], file=paste0(fName,"_",name, "_Q",cutoff, ".txt"), sep="\t", quote = FALSE, row.names = FALSE)

# #now filter based on genotype probabilities
# #implementing this the lazy way now... since I'm looking at things in common between set a and set b, I'm just going to filter based on the first element of each set. Right now I'm using a pretty stringent cutoff, so it's fine
# #use a massive one-liner to split the sample information for one column and 
# apply(subVcf[,setA[1],drop=FALSE],1,function(x) {any(as.numeric(strsplit(strsplit(x,":")[[1]][[8]],",")[[1]])<cutoff)}) -> truthsA
# apply(subVcf[,setB[1],drop=FALSE],1,function(x) {any(as.numeric(strsplit(strsplit(x,":")[[1]][[8]],",")[[1]])<cutoff)}) -> truthsB
# #now further subset the data
# subVcf2 <- subVcf[which(truthsA & truthsB),]




  



