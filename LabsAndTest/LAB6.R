install.packages('seqinr')
install.packages('Biostrings')
load('Biostrings')
load('seqinr')

PUN03.1 <- readDNAStringSet('FCA877.fasta')
PUN03.2 <- readDNAStringSet('FCA883.fasta')
PUN03.3 <- readDNAStringSet('FCA894.fasta')
PUN03.4 <- readDNAStringSet('FCA917.fasta')
PUN03.5 <- readDNAStringSet('FCA924.fasta')

PUNSeq <- c(PUN03.1, PUN03.2, PUN03.3, PUN03.4, PUN03.5)
BiocManager::install('msa')
Library('msa')
Muscle_Alignment <- msa(PUNSeq, type = 'dna', method = 'Muscle')
Muscle_Alignment
print(Muscle_Alignment, show='complete')
Muscle_Alignment<- as.matrix(Muscle_Alignment)
sum(Muscle_Alignment == "-")
aaseq <- readDNAStringSet('FCA877.fasta')
dna1 <- aaseq[[1]]
aaseq1 <- Biostrings::translate(dna1)
install.packages(Phangorn)
Library(Phangorn)
Alignment_phyDat <- msaConvert(Muscle_Alignment, type="phangorn::phyDat")
> write.phyDat(Alignment_phyDat, "alignment.fasta", format = "fasta")
allseq <- toupper(unlist(strsplit(paste(alignment$seq, collapse=""), "")))
allseq <- allseq[allseq %in% c("A","T","G","C")]
gc_total <- sum(allseq == "G" | allseq == "C") / length(allseq)
gc_total * 100

seq_list <- lapply(Muscle_Alignment, function(x) strsplit(x, "")[[1]])
names(seq_list) <- Muscle_Alignment
dna <- as.DNAbin(seq_list)

dist_matrix <- dist.dna(dna, model = "raw")
print(as.matrix(dist_matrix))