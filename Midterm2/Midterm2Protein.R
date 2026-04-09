library(Biostrings)

library(UniprotR)
library(protti)
library(r3dmol)
library(seqinr)
# Some of these may not be needed I'm copying from my lab10
Midterm2Align <- readDNAStringSet('metazoa_alignment.gene.fasta')
# Feel like this is self explanatory
names(Midterm2Align)
# Not really needed but shows names of each sequence in alignment to idenify 'Homo sapiens'
gaps <- letterFrequency(Midterm2Align, '-')
# Identify Homo sapiens
hs_seq <- Midterm2Align[grep("Homo", names(Midterm2Align))]
hs_seq
writeXStringSet(hs_seq, "Homo_sapiens_gene.fasta")
# Pretty much isolates the human sequence from the alignment and writes it to it's own fasta file
#Accession number ncbi X98093 uniprot P54098
accession <- 'P54098'
fetch_uniprot(accession)
AccGoInfo <- GetProteinGOInfo(accession)
PlotGoInfo(AccGoInfo)
# All 4 steps found accession wrote it to a command used uniprot commands to fetch accession from database to allow me to fetch GO terms and then plot.
