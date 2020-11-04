## PastDB scripts

Code associated with the PastDB publication (Martin et al, 2020).

------

* Pipeline to obtain abiotic and stress core AS sets:


* Pipeline to obtain abiotic and stress core GE sets:

  
* Scripts in bin (all perl scripts contain a help option on how to be run):

  - Get_Event_Stats.pl: to calculate general statistics per AS event from any INCLUSION table.
  - Get_PanAS_Events.pl: to define PanAS events from any INCLUSION table.
  - Get_Tissue_Specific_AS.pl: to get tissue-specific AS events from any INCLUSION table.
  - Get_Tissue_Specific_GE.pl: to get genes with tissue-specific expression from any cRPKM/TPM table.
  - Quantify_AS_by_Subsampling.pl: calculate the fraction of genes that are alternatively spliced by event type from an INCLUSION table.
  - Calculate_SS_SCORES_From_PWMs.R: to calculate PWM-based splice site scores.
  - Pipeline_Get_Chain_Aln.sh: bash pipeline to obtain liftOver files.
  - Get_Results_From_Liftover.pl: used to parse the pairwise liftover outputs
  - Get_Results_From_ExOrthist.pl: used to perform the 4-way overlap between core AS sets.


* Files in data folder:

  - General files:
    - AllEvents_for_comparison-Ath.txt.gz (1.3M)
    - Ath.Event-Gene.IDs.txt (9.7M)

  - Config file for tissue-specific analyses:
    - config_for_TS_Ath.txt (551B)

  - Splice sites to calculate SS scores based on PWMs:
    - Annotated_ACCEPTORS-Ath.fasta.gz (1.2M)
    - Annotated_DONORS-Ath.fasta.gz (615K)
    - REFERENCE-ALL_ANNOT-Ath163-3ss.fasta.gz (3.6M)
    - REFERENCE-ALL_ANNOT-Ath163-5ss.fasta.gz (1.8M)
  
  - Lifted events to Brassicacea species by event type:
    - EX-Ath-to-Aal-FILTERED.tab.gz (803K)
    - EX-Ath-to-Aly-FILTERED.tab.gz (920K)
    - EX-Ath-to-Bra-FILTERED.tab.gz (766K)
    - EX-Ath-to-Csa-FILTERED.tab.gz (892K)
    - INT-Ath-Aal-FILTERED.tab.gz (498K)
    - INT-Ath-Aly-FILTERED.tab.gz (1.0M)
    - INT-Ath-Bra-FILTERED.tab.gz (573K)
    - INT-Ath-Csa-FILTERED.tab.gz (941K)
    - ALTA-Ath-to-Aal-FILTERED.tab.gz (459K)
    - ALTA-Ath-to-Aly-FILTERED.tab.gz (648K)
    - ALTA-Ath-to-Bra-FILTERED.tab.gz (412K)
    - ALTA-Ath-to-Csa-FILTERED.tab.gz (580K)
    - ALTD-Ath-to-Aal-FILTERED.tab.gz (245K)
    - ALTD-Ath-to-Aly-FILTERED.tab.gz (356K)
    - ALTD-Ath-to-Bra-FILTERED.tab.gz (216K)
    - ALTD-Ath-to-Csa-FILTERED.tab.gz (315K)

  - Gene and exon orthology clusters:
    - gene_cluster_file-araTha10_ce11_dm6_hg38.gz (286K)
    - EX_clusters-int2b.tab (3.3M)


