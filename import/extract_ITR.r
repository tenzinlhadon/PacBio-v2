#extract ITRs from 3 prime
getwd()
args = commandArgs(trailingOnly=TRUE)
library(GenomicAlignments)
library(rjson)
ref_vec = readLines("ref_heterogeneity.txt")
sample = args[1]
ITR_3_prime = args[2]
ITR_5_prime = args[3]
extract_ITR <- function(ref,three_prime_ITR,five_prime_ITR, length){
    file_path_3_prime = paste(sample,length,  "3_prime.fa", sep = "_")
    file_path_5_prime = paste(sample,length,  "5_prime.fa",sep = "_")
    file_bam1= paste(sample,"_heterogeneity.", length,".sorted.bam",sep ="")
    ITR_length_3_prime = paste(ref,":",three_prime_ITR, sep = "")
    ITR_length_3_prime = gsub('[\"]',"", ITR_length_3_prime)
    ITR_length_5_prime = paste(ref,":",five_prime_ITR, sep = "")
    ITR_length_5_prime = gsub('[\"]',"", ITR_length_5_prime)
    stack_3_prime <- stackStringsFromBam(file_bam1, param=GRanges(ITR_length_3_prime), use.names = T )
    stack_5_prime <- stackStringsFromBam(file_bam1, param=GRanges(ITR_length_5_prime), use.names = T )
    writeXStringSet(stack_3_prime , filepath = file_path_3_prime, format = "fasta", append = TRUE)
    writeXStringSet(stack_5_prime , filepath = file_path_5_prime, format = "fasta", append = TRUE)
}
for(ref in ref_vec){
extract_ITR(ref,ITR_3_prime,ITR_5_prime,"intact")
extract_ITR(ref,ITR_3_prime,ITR_5_prime,"partial")      
}

