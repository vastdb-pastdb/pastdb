#!/usr/bin/sh
# UCSC utilities can be downloaded from: http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/
# G1 is the root "Ath" or "Csa"
# Needs to create symbolic links to the genes (G1.fasta and G2.fasta)
G1=$1
G2=$2
ID_score="80"
 
# split G2 into individual chromosome/scaffold (G2.$i.fa) and save them into ./G2_chr direcotry
mkdir -p CHR_FA
perl ../SplitByChr_forSh.pl $G2

n=$(ls ./CHR_FA/$G2*fa | wc -l)
 
# Split the G2 chromosomes/scaffolds into 3K chunks and make lift files
mkdir -p lift
mkdir -p split
for i in $(seq 1 $n)
do
  faSplit size ./CHR_FA/$G2.$i.fa  3000 ./split/$G2.$i.split  -lift=./lift/$G2.$i.lft -oneFile
done
 
# run blat (def was -minIdentity=95 -minScore=100)
mkdir -p psl
for i in $(seq 1 $n)
do
  blat $G1.fasta ./split/$G2.$i.split.fa -t=dna -q=dna -tileSize=12 -fastMap -minIdentity=$ID_score -noHead -minScore=50  ./psl/$G2.$i.psl
done
 
# Change coordinates of .psl files to parent coordinate system
mkdir -p liftup
for i in $(seq 1 $n)
do
    liftUp -pslQ ./liftup/$G2.$i.liftup.psl ./lift/$G2.$i.lft warn ./psl/$G2.$i.psl
done
 
# Make chain files
mkdir -p chain_raw
for i in $(seq 1 $n)
do
    axtChain -linearGap=medium -faQ -faT -psl ./liftup/$G2.$i.liftup.psl ./$G1.fasta ./$G2.fasta ./chain_raw/$i.chain
done
 
# Merge and sort chain files
chainMergeSort ./chain_raw/*.chain | chainSplit chain_split stdin
 
faSize $G1.fasta  -detailed >$G1.chr_length.txt
faSize $G2.fasta  -detailed >$G2.chr_length.txt
 
# Make alignment nets from chain files
mkdir -p net
for i in  ./chain_split/*.chain
do
 tag=${i/\.\/chain_split\//}
 chainNet $i ./$G1.chr_length.txt ./$G2.chr_length.txt ./net/$tag.net /dev/null
done
 
# Create liftOver chain file
mkdir -p over
for i in ./chain_split/*.chain
do
  tag=${i/\.\/chain_split\//}
  netChainSubset  ./net/$tag.net $i ./over/$tag.chain
done

To='To'
 
cat ./over/*.chain >$G1$To$G2.over.chain
