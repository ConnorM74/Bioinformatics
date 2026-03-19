setwd()
install.packages('UniprotR')
install.packages('protti')
install.packages('r3dmol')
BiocManager::install('GenomicAlignments')
library(UniprotR)
library(protti)
library(r3dmol)
library(seqinr)
library(Biostrings)
library(GenomicAlignments)
ProSeq <- readAAStringSet('Hsapiens6.fasta')
Biostrings::writeXStringSet(ProSeq)
CsvProSeq <- read.csv('ProSeqLAB10.txt')
ChProSeq <- c('P0A799', 'P08839')
ProSeqGoInfo <- GetProteinGOInfo(ChProSeq)
PlotGoInfo(ProSeqGoInfo)
# Below is code sourced from github as stated in lab to save graphs.
PlotGOAll(GOObj = ProSeqGoInfo,
  Top = 10,
  directorypath = getwd(),
  width = 8, 
  height = 5)
GetPathology_Biotech(ProSeqGoInfo) 
#This command just isn't working whatsoever I used Ai to give me this command 
accessions <- rownames(ProSeqGoInfo)  # if your data frame rows are accessions

Pathology_list <- lapply(accessions, function(x) {tryCatch(GetPathology_Biotech(x),error = function(e) NULL)})
PathologyInfo <- do.call(rbind, Pathology_list[!sapply(Pathology_list, is.null)])
#Just an alternative way to get pathology info desperate attempt at this point.

PathologyInfo
Get.diseases(ProSeqGoInfo)
#Continuing to have issues with every single command switching to doing this manually
fetch_uniprot('P0A799')
fetch_uniprot('P08839')
fetch_pdb('1ZMR')
fetch_pdb('2HWG')
fetch_alphafold_prediction('1ZMR')
fetch_alphafold_prediction('2HWG')
