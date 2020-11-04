if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("Biostrings")
BiocManager::install("XVector")

library(Biostrings)
library(XVector)

### For 3 ss
# load fasta with reference SS
fasr <- readDNAStringSet("Annotated_ACCEPTORS-Ath.fasta", format= "fasta")
# load fasta with all SS
fasi <- readDNAStringSet("REFERENCE-ALL_ANNOT-Ath163-3ss.fasta", format= "fasta")
# create PWM from ref fasta
pfmi <- consensusMatrix(toupper(fasr), as.prob=F) # create frequency matrix (counts)
pwmi <- PWM(pfmi) # PWM from freq matrix
pwmi <- pwmi[c("A","C","G","T"),] # clean
# score all SS relative to PWM obtained from ref fasta
scoi <- sapply(fasi, function(x) { PWMscoreStartingAt(pwmi,x,starting.at = 1) })
# save
scod <- as.data.frame(scoi)
write.table(scod,col.names=F,file="REFERENCE-ALL_ANNOT-Ath163-3ss.PWM.scores", quote=F, sep="\t")

### For 5 ss
# load fasta with reference SS
fasr <- readDNAStringSet("Annotated_DONORS-Ath.fasta", format= "fasta")
# load fasta with all SS
fasi <- readDNAStringSet("REFERENCE-ALL_ANNOT-Ath163-5ss.fasta", format= "fasta")
# create PWM from ref fasta
pfmi <- consensusMatrix(toupper(fasr), as.prob=F) # create frequency matrix (counts)
pwmi <- PWM(pfmi) # PWM from freq matrix
pwmi <- pwmi[c("A","C","G","T"),] # clean
# score all SS relative to PWM obtained from ref fasta
scoi <- sapply(fasi, function(x) { PWMscoreStartingAt(pwmi,x,starting.at = 1) })
# save
scod <- as.data.frame(scoi)
write.table(scod,col.names=F,file="REFERENCE-ALL_ANNOT-Ath163-5ss.PWM.scores", quote=F, sep="\t")
