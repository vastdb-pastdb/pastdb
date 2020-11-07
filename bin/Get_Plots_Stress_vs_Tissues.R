### Script used to plot Figure 5c
rm (list=ls())
Table=read.table("Stress_vs_Tissues-input_table.tab",header=T, row.names=NULL, sep="\t")
# Input format:
# 1: EventID
# 2: PSI range across tissues
# 3: Maximum dPSI among abiotic stress experiments
# 4: Maximum dPSI among biotic stress experiments
# 5: N of tissues with coverage
# 6: N of abiotic stress experiments with coverage
# 7: N of biotic stress experiments with coverage
# 8: Species

# It calculates the global PSI variation, and the relative contributions
Table$GlobalPSI <- Table$Range_tis + Table$Max_dPSI_Abiotic + Table$Max_dPSI_biotic
Table$T_prop <- (Table$Range_tis/Table$GlobalPSI)*100
Table$A_prop <- (Table$Max_dPSI_Abiotic/Table$GlobalPSI)*100
Table$B_prop <- (Table$Max_dPSI_biotic/Table$GlobalPSI)*100

# Cut-offs
i <- 10 # min sum of the three values (tissue range + max stress dPSIs)
j <- 4  # min number of tissues/experiments with coverage in each type

# Data loading
Ath_T <- na.omit(as.vector(Table[Table$Species=="Ath" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,10]))
Ath_A <- na.omit(as.vector(Table[Table$Species=="Ath" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,11]))
Ath_B <- na.omit(as.vector(Table[Table$Species=="Ath" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))
Cel_T <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,10]))
Cel_A <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,11]))
Cel_B <- na.omit(as.vector(Table[Table$Species=="Cel" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))
Dme_T <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,10]))
Dme_A <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,11]))
Dme_B <- na.omit(as.vector(Table[Table$Species=="Dme" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))
Hsa_T <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,10]))
Hsa_A <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,11]))
Hsa_B <- na.omit(as.vector(Table[Table$Species=="Hsa" & Table$GlobalPSI>i & Table$N_tis>=j & Table$N_abiotic>=j & Table$N_biotic>=j,12]))

#HEXBIN PLOT
library(hexbin)
library(pals)

binAth_A <- (hexbin(Ath_A,Ath_T,xbins = 20))
binAth_B <-  (hexbin(Ath_B,Ath_T,xbins = 20))
binCel_A <- (hexbin(Cel_A,Cel_T,xbins = 20))
binCel_B <- (hexbin(Cel_B,Cel_T,xbins = 20))
binDme_A  <- (hexbin(Dme_A,Dme_T,xbins = 20))
binDme_B <- (hexbin(Dme_B,Dme_T,xbins = 20))
binHsa_A <- (hexbin(Hsa_A,Hsa_T,xbins = 20))
binHsa_B <-(hexbin(Hsa_B,Hsa_T,xbins = 20))

plot (binAth_A, main="" , colramp=coolwarm)
plot (binAth_B, main="" , colramp=coolwarm)
plot (binCel_A, main="" , colramp=coolwarm)
plot (binCel_B, main="" , colramp=coolwarm)
plot (binDme_A, main="" , colramp=coolwarm)
plot (binDme_B, main="" , colramp=coolwarm)
plot (binHsa_A, main="" , colramp=coolwarm)
plot (binHsa_B, main="" , colramp=coolwarm)


