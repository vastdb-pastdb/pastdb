# PastDB scripts

Code associated with the PastDB publication (Martin et al, 2020) 

* Pipeline to obtain abiotic and stress core AS sets:
  - 


* Pipeline to obtain abiotic and stress core GE sets:

  

* Scripts in bin:
  - Get_Event_Stats.pl: to calculate general statistics per AS event from an INCLUSION table.
  - Get_PanAS_Events.pl: to define PanAS events.
  - Get_Tissue_Specific_AS.pl: to get tissue-specific AS events.
  - Get_Tissue_Specific_GE.pl: to get genes with tissue-specific expression.
  - Quantify_AS_by_Subsampling.pl: calculate the fraction of genes that are alternatively spliced by event type.
  - Calculate_SS_SCORES_From_PWMs.R: to calculate PWM-based splice site scores.
  - Pipeline_Get_Chain_Aln.sh: bash pipeline to obtain liftOver files.
  - Get_Results_From_Liftover.pl: used to parse the pairwise liftover outputs
  - Get_Results_From_ExOrthist.pl: used to perform the 4-way overlap between core AS sets.


* Files in data folder:
  - General files:
1.3M	AllEvents_for_comparison-Ath.txt.gz
9.7M	Ath.Event-Gene.IDs.txt

  - Config file for tissue-specific analyses:
551B	config_for_TS_Ath.txt

  - Splice sites to calculate SS scores based on PWMs:
1.2M	Annotated_ACCEPTORS-Ath.fasta.gz
615K	Annotated_DONORS-Ath.fasta.gz
3.6M	REFERENCE-ALL_ANNOT-Ath163-3ss.fasta.gz
1.8M	REFERENCE-ALL_ANNOT-Ath163-5ss.fasta.gz
  
  - Lifted events to Brassicacea species by event type:
803K	EX-Ath-to-Aal-FILTERED.tab.gz
920K	EX-Ath-to-Aly-FILTERED.tab.gz
766K	EX-Ath-to-Bra-FILTERED.tab.gz
892K	EX-Ath-to-Csa-FILTERED.tab.gz
498K	INT-Ath-Aal-FILTERED.tab.gz
1.0M	INT-Ath-Aly-FILTERED.tab.gz
573K	INT-Ath-Bra-FILTERED.tab.gz
941K	INT-Ath-Csa-FILTERED.tab.gz
459K	ALTA-Ath-to-Aal-FILTERED.tab.gz
648K	ALTA-Ath-to-Aly-FILTERED.tab.gz
412K	ALTA-Ath-to-Bra-FILTERED.tab.gz
580K	ALTA-Ath-to-Csa-FILTERED.tab.gz
245K	ALTD-Ath-to-Aal-FILTERED.tab.gz
356K	ALTD-Ath-to-Aly-FILTERED.tab.gz
216K	ALTD-Ath-to-Bra-FILTERED.tab.gz
315K	ALTD-Ath-to-Csa-FILTERED.tab.gz

  - Gene and exon orthology clusters:
286K	gene_cluster_file-araTha10_ce11_dm6_hg38.gz
3.3M	EX_clusters-int2b.tab
