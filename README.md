## PastDB scripts and associated data

Code and data associated with the PastDB web and publication (Martin et al, 2020). For any further enquires, please feel free to contact Manuel Irimia (mirimia@gmail.com) and/or Guiomar Martin (guiomarm@igc.gulbenkian.pt)

------

* Pipeline to obtain abiotic and stress core AS sets:


* Pipeline to obtain abiotic and stress core GE sets:

  
* Scripts in bin (all perl scripts contain a help option on how to be run):

  - `Get_Event_Stats.pl`: to calculate general statistics per AS event from any INCLUSION table.
  - `Get_PanAS_Events.pl`: to define PanAS events from any INCLUSION table.
  - `Get_Tissue_Specific_AS.pl`: to get tissue-specific AS events from any INCLUSION table.
  - `Get_Tissue_Specific_GE.pl`: to get genes with tissue-specific expression from any cRPKM/TPM table.
  - `Quantify_AS_by_Subsampling.pl`: calculate the fraction of genes that are alternatively spliced by event type from an INCLUSION table.
  - `Calculate_SS_SCORES_From_PWMs.R`: to calculate PWM-based splice site scores.
  - `Pipeline_Get_Chain_Aln.sh`: bash pipeline to obtain liftOver files.
  - `Get_Results_From_Liftover.pl`: used to parse the pairwise liftover outputs
  - `Get_Results_From_ExOrthist.pl`: used to perform the 4-way overlap between core AS sets.


* Files from PastDB: the main data files used for the analyses are available for [download in PastDB](http://pastdb.crg.eu/wiki/Downloads), and are also copied here:

  - AS events:
    - [EVENTS table](http://vastdb.crg.eu/downloads/araTha10/EVENT_INFO-araTha10.tab.gz): Information about AS event coordinates and sequences. TAIR10 asssembly.
    - [MAIN PSI table](http://vastdb.crg.eu/downloads/araTha10/PSI_TABLE-araTha10.tab.gz): Inclusion patterns of AS events across tissues, cell types and developmental stages (main PSI plot).
    - [ABIOTIC PSI table](http://vastdb.crg.eu/downloads/araTha10/PSI_TABLE-araTha10-40-ABIOTIC-v251.tab.gz): Inclusion patterns of AS events in ABIOTIC stress experiments (special dataset).
    - [BIOTIC PSI table](http://vastdb.crg.eu/downloads/araTha10/PSI_TABLE-araTha10-18-BIOTIC-v251.tab.gz): Inclusion patterns of AS events in BIOTIC stress experiments (special dataset).
    - [LIGHT PSI table](http://vastdb.crg.eu/downloads/araTha10/PSI_TABLE-araTha10-21-LIGHT-v251.tab.gz): Inclusion patterns of AS events in LIGHT experiments (special dataset).
    - [SPL_FACTORS PSI table](http://vastdb.crg.eu/downloads/araTha10/PSI_TABLE-araTha10-33-SPL_FACTORS-v251.tab.gz): Inclusion patterns of AS events upon SPLICING FACTOR disruption (special dataset).
    
  - Event features:
    - [SPLICE SITES table](http://vastdb.crg.eu/downloads/araTha10/SPLICE_SITE_SCORES-araTha10.tab.gz): Sequences and strength scores of 5' and 3' splice sites of alternative exons.
    - [PCR VALIDATION table](http://vastdb.crg.eu/downloads/araTha10/PCR_PRIMERS-araTha10.tab.gz): Suggested primer sequences and expected band lengths for validation of AS events by RT-PCR.
    
  - Protein impact:
    - [PROTEIN IMPACT table](http://vastdb.crg.eu/downloads/araTha10/PROT_IMPACT-araTha10-v3.tab.gz): Effect of the AS event in the open reading frame of the transcript. Version v3.
    - [PROTEIN ISOFORMS](http://vastdb.crg.eu/downloads/araTha10/PROT_ISOFORMS-araTha10.tab.gz): Mappings of events to ProteinIDs.
    - [DOMAINS table (PFAM)](http://vastdb.crg.eu/downloads/araTha10/PROT_PFAM-araTha10.tab.gz): Mappings to Pfam domains
    - [DOMAINS table (PROSITE)](http://vastdb.crg.eu/downloads/araTha10/PROT_PROSITE-araTha10.tab.gz): Mappings to PROSITE domains.
    - [PROTEIN DISORDERED REGIONS table](http://vastdb.crg.eu/downloads/araTha10/PROT_DISORDER-araTha10.tab.gz): Intrinsic disorder rates for A, C1 and C2 exons, using disopred3.
    
  - Genes:
    - [GENES table](http://vastdb.crg.eu/downloads/araTha10/GENE_INFO-araTha10.tab.gz): Information about gene names, descriptions, genomic coordinates and biotypes.
    - [MAIN EXPRESSION table](http://vastdb.crg.eu/downloads/araTha10/EXPRESSION_TABLE-araTha10.tab.gz): Gene expression across tissues, cell types and developmental stages. Measured in cRPKM and in raw reads (main GE plot).
    - [ABIOTIC EXPRESSION table](http://vastdb.crg.eu/downloads/araTha10/cRPKM-araTha10-40-ABIOTIC-NORM.tab.gz): Gene expression in ABIOTIC stress experiments (special dataset).
    - [BIOTIC EXPRESSION table](http://vastdb.crg.eu/downloads/araTha10/cRPKM-araTha10-18-BIOTIC-NORM.tab.gz): Gene expression in BIOTIC stress experiments (special dataset).
    - [LIGHT EXPRESSION table](http://vastdb.crg.eu/downloads/araTha10/cRPKM-araTha10-21-LIGHT-NORM.tab.gz): Gene expression in LIGHT experiments (special dataset).
    - [SPL_FACTORS EXPRESSION table](http://vastdb.crg.eu/downloads/araTha10/cRPKM-araTha10-33-SPL_FACTORS-NORM.tab.gz): Gene expression upon SPLICING FACTOR disruption (special dataset).
    - [GENE-EVENTS table](http://vastdb.crg.eu/downloads/araTha10/EVENTID_to_GENEID-araTha10.tab.gz): Table relating genes to AS events.

  - Samples:
    - [SAMPLE_INFO table](http://vastdb.crg.eu/downloads/araTha10/SAMPLE_INFO-araTha10.tab.gz): SRA identifiers and other information related to RNA-seq data used in this database.
  


* Files in data/ folder:

  - General files:
    - AllEvents_for_comparison-Ath.txt.gz (1.3M)
    - Ath.Event-Gene.IDs.txt (9.7M)

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


