#install.packages('seqinr')
#install.packages('Biostrings')
#Only install packages if you have not already.
library('Biostrings')
library('seqinr')
#Library activates installed packages
MSeqs <- readDNAStringSet('MidtermSeqs.fasta')
MSeqs
#Allows us to view fasta file as a DNA string set.


#BiocManager::install('msa')
Library('msa')
MidMSA <- msa(MSeqs, type = 'dna', method = 'Muscle')
MidMSA
#Allows us to create a MSA on our fasta file sequences.
#Also generates a "Consensus sequence.
MidCon <- msaConsensusSequence(MidMSA)
#Saves consesus sequence
AlnDNA <- as(MidMSA, 'DNAStringSet')
MidPid <- 1 - (as.matrix(dist_mat) / width(AlnDNA)[1])
mean(MidPid)
#This changes our MSA's class to DNAStringSet
#After we change it to a DNAStringSet we take the distance matrix
# and divide it by the width of the Alignment to get our Pid
# Pid [1] 0.998528 = 99.8%
MidConVec <- strsplit(MidCon, "")[[1]]  
gc <- sum(MidConVec %in% c("G", "C"))  
gc_content <- gc / length(MidConVec)  
gcPercent <- gc_content*100
#This turns the sequence into readable individual bases.
#I then summed up all the G's and C's, divded them by the length.
#This gives gc in fractional for so x100 and got 51.6%
MatDNA <- as.matrix(AlnDNA)
SeqDiff <- MatDNA != matrix(rep(MidConVec, nrow(MatDNA)),
nrow = nrow(MatDNA), byrow = TRUE)
#Sorts by row to determine if each base pair of 
#each sequence matches the consesus sequence, TRUE= DIFFERENCE.
SeqDiffWhere <- apply(SeqDiff, 1, function(x) which(x))  
SeqDiffWhere
#This command tells you where to look to find your differing bases.
#There are multiple mutations. Some seem to be potential deletions denoted from gaps, and substitutions. 
#Insertions also possible. Transition implied at position 47, 134, and 623. Transversions implied at position 145, 152, 586, and
# QUESTION 6 ANSWER: USING NCBI BLAST
# ACCESSION NUMBER: LC121775.1 
# GENE: Select seq LC121775.1	Homo sapiens hbb gene for beta globin, partial cds, note: HbLimassol Cd8(AAG>AAC)

# Question 7 and on. Hsapiens 6 is the most different
MostDiff <- AlnDNA["Homo_sapiens_6"]
MostDiff1 <- DNAString(gsub("-", "", as.character(MostDiff)))
#Strangely a random gap was found in the sequence but not in the fasta file.
#This implies a deletion at some point if the overall sequence is shorter and supplemented by gaps.
ProAln <- Biostrings::translate(MostDiff1)
#Simple Nucleotide -> Protein translation.
ProAln1 <- AAStringSet(ProAln)
names(ProAln1) <- 'Homo_sapiens_6protein'
writeXStringSet(ProAln1, filepath = "Hsapiens6.fasta")
#In order to create the fasta file i had to convert ProAln from AAstring -> AAstringset.

#Question 8: ACCESSION NUMBER: KAI2558340.1

#Question 9: This gene is responsible for diseases such as B-thalassemia, sickle cell, and other unique mutations.
#The individual identified in question 7, Hsapiens6, likely has B-thalassemia.
#We know this because this disease usually produces non-functional B subunits and
#it's observed that in the case of Hsapiens6, the start codon is AAT which is not a universally accepted start codon.
#Likely this leads to a non-functional protein associated with B-Thalassemia.
