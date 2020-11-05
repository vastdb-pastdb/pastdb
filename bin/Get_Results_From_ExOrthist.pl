#!/usr/bin/perl
# Original script: MatchExonOrths.pl
use strict;
use warnings;

# 1) it processes orthologous gene clusters
open (CLUSTERS, "gene_cluster_file-araTha10_ce11_dm6_hg38") || die "Gene clusters\n";
my %g_cl; my %cl_g;
while (<CLUSTERS>){
    chomp;
    my @t=split(/\t/);
    $g_cl{$t[1]}{$t[2]}=$t[0];
    push(@{$cl_g{$t[0]}{$t[1]}},$t[2]);
}
close CLUSTERS;

### 2) it processes orthologous exon clusters 
open (EX_CLUSTERS, "EX_clusters-int2b.tab") || die "Exon clusters\n";
my %ex_cl; my %cl_ex; 
my %all_ex_cl_per_g_cl;
my %co_ex_cl;
<EX_CLUSTERS>;
while (<EX_CLUSTERS>){
    chomp;
    my @t=split(/\t/);
    $ex_cl{$t[3]}{$t[2]}=$t[0]; # needs species (t[3]), since the exons are only coordinates, not IDs
    push(@{$cl_ex{$t[0]}{$t[3]}},$t[2]);
    my ($g_cl)=$t[0]=~/(.+?)\./;
    push(@{$all_ex_cl_per_g_cl{$g_cl}},$t[0]);

    ### to match to VTS events later on; it stores the 3' and the 5'ss
    my ($chr,$i,$f)=$t[2]=~/(.+?)\:(.+?)\-(.+)/;
    my $coA="$chr:$i";
    my $coB="$chr:$f";
    $co_ex_cl{$t[3]}{$coA}=$t[0];
    $co_ex_cl{$t[3]}{$coB}=$t[0];
}
close EX_CLUSTERS;


my @SPECIES=("araTha10","ce11","dm6","hg38");
my %is_EX_5UTR;
my %is_G_cl_core;
my %is_G_cl_core_5UTR;
my %is_EX_cl_core;
my %ev_g;
my %cores_ev;
my %has_cluster; my %has_cluster1; my %has_cluster2;

foreach my $sp (@SPECIES){
    ### 3) Loads ONTO information. Mainly to analyze 5' UTR events separately
    my $onto_file = "ONTOs/$sp"."_ONTO_all-v3.tab";
    open (ONTO, $onto_file) || die "Cannot open $onto_file\n";
    while (<ONTO>){
	chomp;
	my @t=split(/\t/);
	$is_EX_5UTR{$t[1]}=1 if $t[7]=~/UTR_5/; # simply stores which events are 5'UTRs
    }
    close ONTO;

    ### 4) Loads core set information for each species and process the events
    open (CORES, "Core_$sp-ALL.tab") || die "Cannot open Core_$sp-ALL.tab\n";
    while (<CORES>){
	chomp;
	my @t=split(/\t/);
	my $ev = $t[1];
	my $g=$t[6];
	my $core=$t[7];
	my ($type)=$ev =~ /.{3}(.+?)\d/;
	
	next if $core eq "AS_NC";

	### 4.1) Gene level match
	if ($g_cl{$sp}{$g}){ # checks if there is a cluster for that gene
	    $is_G_cl_core{$g_cl{$sp}{$g}}{$sp}{$core}.="$ev,"; # adds all events for all types for that gene cluster, per CORE.
	    $is_G_cl_core_5UTR{$g_cl{$sp}{$g}}{$sp}{$core}.="$ev," if $is_EX_5UTR{$ev}; # i.e. if there is at least one gene that has a 5 UTR event, it will exist
	}

	### 4.2) Exon level match
	$ev_g{$ev}=$g;
	push(@{$cores_ev{$core}{$sp}},$ev); # keeps all events for each CORE type and SPECIES for later
	
	### 4.2.1) It analyzes each event type separately, to see if any matches to an exon in an orthologous exon cluster
	###        It does it by trying to match either of both splice sites, depending on the event type.
	my $cl_match=""; # EX, ALTA, ALTD
	my $cl_match1=""; my $cl_match2=""; # INT

	if ($type eq "EX"){ # matches both splice sites
	    my ($chr,$I,$F)=$t[4]=~/(.+?)\:.+?\,(.+?)\-(.+?)\,/;
	    my @I=split(/\+/,$I);
	    my @F=split(/\+/,$F);
	    foreach my $temp (@I){
		my $co="$chr:$temp";
		$cl_match = $co_ex_cl{$sp}{$co} if !$cl_match && $co_ex_cl{$sp}{$co};
	    }
	    foreach my $temp (@F){
		my $co="$chr:$temp";
		$cl_match = $co_ex_cl{$sp}{$co} if !$cl_match && $co_ex_cl{$sp}{$co};
	    }
	    if ($cl_match){ # EX
		$has_cluster{$ev}=$cl_match;
		$is_EX_cl_core{$cl_match}{UP}{$sp}{$core}.="$ev,"; # all kinds of AS types matching that exon cluster
		$is_EX_cl_core{$cl_match}{DOWN}{$sp}{$core}.="$ev,"; # all kinds of AS types matching that exon cluster
	    }
	}
	elsif ($type eq "INT"){ # it matches the 5'ss of the upstream exon and the 3'ss of the downstream
	    my ($chr,$i1,$f1,$i2,$f2,$str)=$t[4]=~/(.+?)\:(.+?)\-(.+?)\=(.+?)\-(.+?)\:(.)/;
	    if ($str eq "+"){
		my $coA="$chr:$f1"; my $coB="$chr:$i2";
		$cl_match1 = $co_ex_cl{$sp}{$coA} if !$cl_match1 && $co_ex_cl{$sp}{$coA};
		$cl_match2 = $co_ex_cl{$sp}{$coB} if !$cl_match2 && $co_ex_cl{$sp}{$coB};
	    }
	    elsif ($str eq "-"){
		my $coA="$chr:$i1"; my $coB="$chr:$f2";
		$cl_match1 = $co_ex_cl{$sp}{$coA} if !$cl_match1 && $co_ex_cl{$sp}{$coA};
		$cl_match2 = $co_ex_cl{$sp}{$coB} if !$cl_match2 && $co_ex_cl{$sp}{$coB};    
	    }
	    if ($cl_match1){ # INT (upstream exon is matched in its 5'ss [DOWN]
		$has_cluster1{$ev}=$cl_match1;
		$is_EX_cl_core{$cl_match1}{DOWN}{$sp}{$core}.="$ev,"; # all kinds of AS types matching that exon cluster
	    }
	    if ($cl_match2){ # INT (downstream exon is matched in its 3'ss [UP]
		$has_cluster2{$ev}=$cl_match2;
		$is_EX_cl_core{$cl_match2}{UP}{$sp}{$core}.="$ev,"; # all kinds of AS types matching that exon cluster
	    }
	}
	elsif ($type eq "ALTA"){ # only tries to match the 3'ss of the exon in the clusters (UP)
	    my ($chr,$I,$F)=$t[4]=~/(.+?)\:.*?\,(.*)\-(.*)/;
	    my $str="";
	    $str="+" if $I=~/\+/;
	    $str="-" if $F=~/\+/;
	    print "ISSUE: no $str for $ev\n" if !$str;

	    my @I=split(/\+/,$I);
	    my @F=split(/\+/,$F);
	    my @CO;
	    @CO=@I if $str eq "+";
	    @CO=@F if $str eq "-";
	    
	    foreach my $temp (@CO){
		my $co="$chr:$temp";
		$cl_match = $co_ex_cl{$sp}{$co} if !$cl_match && $co_ex_cl{$sp}{$co};
	    }
	    if ($cl_match){ # ALTA
		$has_cluster{$ev}=$cl_match;
		$is_EX_cl_core{$cl_match}{UP}{$sp}{$core}.="$ev,"; # all kinds of AS types matching that exon cluster
	    }
	}
	elsif ($type eq "ALTD"){ # only tries to match the 5'ss of the exon in the clusters (DOWN) 
	    my ($chr,$I,$F)=$t[4]=~/(.+?)\:(.*?)\-(.*?)\,/;
	    my $str="";
	    $str="+" if $F=~/\+/;
	    $str="-" if $I=~/\+/;
	    print "ISSUE: no $str for $ev\n" if !$str;

	    my @I=split(/\+/,$I);
	    my @F=split(/\+/,$F);
	    my @CO;
	    @CO=@F if $str eq "+";
	    @CO=@I if $str eq "-";
	    
	    foreach my $temp (@CO){
		my $co="$chr:$temp";
		$cl_match = $co_ex_cl{$sp}{$co} if !$cl_match && $co_ex_cl{$sp}{$co};
	    }
	    if ($cl_match){ # ALTD
		$has_cluster{$ev}=$cl_match;
		$is_EX_cl_core{$cl_match}{DOWN}{$sp}{$core}.="$ev,"; # all kinds of AS types matching that exon cluster
	    }
	}
    }
    close CORES;
}

my @CORES=("ABI","BIO","TIS");
### Main Counters
my %tally_G_orth; # exons with gene orth in sp2
my %tally_G_orth_core; # exons with gene orth in sp2 that also have a core event (any type).
my %tally_EX_orth; # exons with an exon orth in sp2 where event falls 
my %tally_EX_orth_core; # exons with an exon orth in sp2 that also have a core event in the equivalent SS of the orth exon

my %tally_EX_5UTR_all; # number of 5'UTR events in Sp1
my %tally_EX_5UTR; # number of 5'UTR events in Sp1 among gene orth with core events in both species
my %tally_EX_5UTR_core; # number of 5'UTR events in Sp1 for which a gene orth of sp2 also has a 5'UTR core event

my %combined_G_orth_core; # per event to do the 4-way
my %combined_EX_orth_core; # same as tally_EX_orth_core + tally_EX_5UTR_core but 4-way

### 5) It goes event by event and species and it counts the number of matches at different levels + prints them out
foreach my $core (@CORES){
    foreach my $sp1 (@SPECIES){
	my @events=@{$cores_ev{$core}{$sp1}};
	
	foreach my $ev (@events){ # One entry per exon, but the same gene could be multiple times.
	    my ($type)= $ev =~ /.{3}(.+?)\d/;
	    my $g1=$ev_g{$ev};
	    $tally_EX_5UTR_all{$core}{$sp1}++ if $is_EX_5UTR{$ev}; # all events in 5'UTRs in the core sets
	    foreach my $sp2 (@SPECIES){
		if ($sp1 ne $sp2){
		    if ($g_cl{$sp1}{$g1} && ($#{$cl_g{$g_cl{$sp1}{$g1}}{$sp2}}>=0)){ # && $#{$cl_g{$g_cl{$sp1}{$g1}}{$sp2}}<6)){ # between 1 to 6 orths
			$tally_G_orth{$core}{$sp1}{$sp2}++; # the gene with the event has an ortholog (in principle, independently of how big the cluster is)

			if ($is_G_cl_core{$g_cl{$sp1}{$g1}}{$sp2}{$core}){ 
			    $tally_G_orth_core{$core}{$sp1}{$sp2}++; # the gene ortholog also has a core event			    
			    $combined_G_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram

			    if ($type ne "INT"){
				if ($has_cluster{$ev} && ($#{$cl_ex{$has_cluster{$ev}}{$sp2}}>=0)){ # i.e. it's in a cluster for the species 2
				    $tally_EX_orth{$core}{$sp1}{$sp2}++; # the event is in an exon cluster

				    ### I.e. the actual exon cluster is core in BOTH species in the same SS
				    if ($type eq "EX"){
					if ($is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core} || $is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}){
					    $tally_EX_orth_core{$core}{$sp1}{$sp2}++;
					    $combined_EX_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram 
					    
					    # just printing the positive hits
					    if ($is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core} && $is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}){
						$is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core}=~s/\,$//;
						$is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}=~s/\,$//;
						print "$core\t$sp1\t$sp2\t$has_cluster{$ev}\t$ev\t$is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core}\t$is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}\n";
					    }
					    elsif ($is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core} && !$is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}){
						$is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core}=~s/\,$//;
						print "$core\t$sp1\t$sp2\t$has_cluster{$ev}\t$ev\t$is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core}\tNA\n";
					    }
					    elsif (!$is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core} && $is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}){
						$is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}=~s/\,$//;
						print "$core\t$sp1\t$sp2\t$has_cluster{$ev}\t$ev\tNA\t$is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}\n";
					    }
					}
				    }
				    elsif ($type eq "ALTA"){
					if ($is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core}){
					    $tally_EX_orth_core{$core}{$sp1}{$sp2}++;
					    $combined_EX_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram 

					    # just printing the positive hits
					    $is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core}=~s/\,$//;
					    print "$core\t$sp1\t$sp2\t$has_cluster{$ev}\t$ev\t$is_EX_cl_core{$has_cluster{$ev}}{UP}{$sp2}{$core}\n";
					}
				    }
				    elsif ($type eq "ALTD"){
					if ($is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}){
					    $tally_EX_orth_core{$core}{$sp1}{$sp2}++;
					    $combined_EX_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram 
					    
					    # just printing the positive hits
					    $is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}=~s/\,$//;
					    print "$core\t$sp1\t$sp2\t$has_cluster{$ev}\t$ev\t$is_EX_cl_core{$has_cluster{$ev}}{DOWN}{$sp2}{$core}\n";
					}
				    }
				    else { print "Issues with $type for $ev (unexpected type)\n";}
				}
				elsif ($is_EX_5UTR{$ev}){ # special case for UTRs
				    $tally_EX_5UTR{$core}{$sp1}++;
				    if ($is_G_cl_core_5UTR{$g_cl{$sp1}{$g1}}{$sp2}{$core}){ # both are 5' UTR and core
					$tally_EX_5UTR_core{$core}{$sp1}{$sp2}++;
					$combined_EX_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram 

					# just printing the positive hits
					$is_G_cl_core_5UTR{$g_cl{$sp1}{$g1}}{$sp2}{$core}=~s/\,$//;
					print "$core\t$sp1\t$sp2\tUTR5\t$ev\t$is_G_cl_core_5UTR{$g_cl{$sp1}{$g1}}{$sp2}{$core}\n";
				    }
				}
			    }
			    elsif ($type eq "INT"){
				# i.e. either of the exons is in a cluster
				if (($has_cluster1{$ev} && ($#{$cl_ex{$has_cluster1{$ev}}{$sp2}}>=0)) || 
				    ($has_cluster2{$ev} && ($#{$cl_ex{$has_cluster2{$ev}}{$sp2}}>=0))){
				    $tally_EX_orth{$core}{$sp1}{$sp2}++; # the event is in an exon cluster

				    ### I.e. the actual exon cluster is core in BOTH species in the same SS
				    if ($has_cluster1{$ev} && ($is_EX_cl_core{$has_cluster1{$ev}}{DOWN}{$sp2}{$core})){
					$tally_EX_orth_core{$core}{$sp1}{$sp2}++;
					$combined_EX_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram 

					# just printing the positive hits
					$is_EX_cl_core{$has_cluster1{$ev}}{DOWN}{$sp2}{$core}=~s/\,$//;
					print "$core\t$sp1\t$sp2\t$has_cluster1{$ev}\t$ev\t$is_EX_cl_core{$has_cluster1{$ev}}{DOWN}{$sp2}{$core}\n";
				    }
				    elsif ($has_cluster2{$ev} && ($is_EX_cl_core{$has_cluster2{$ev}}{UP}{$sp2}{$core})){ # it does not multi-count if both hit a core ortholog
					$tally_EX_orth_core{$core}{$sp1}{$sp2}++;
					$combined_EX_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram 

					# just printing the positive hits
					$is_EX_cl_core{$has_cluster2{$ev}}{UP}{$sp2}{$core}=~s/\,$//;
					print "$core\t$sp1\t$sp2\t$has_cluster2{$ev}\t$ev\t$is_EX_cl_core{$has_cluster2{$ev}}{UP}{$sp2}{$core}\n";
				    }
				}
				elsif ($is_EX_5UTR{$ev}){ # special case for UTRs
				    $tally_EX_5UTR{$core}{$sp1}++;
				    
				    if ($is_G_cl_core_5UTR{$g_cl{$sp1}{$g1}}{$sp2}{$core}){ # both are 5' UTR and core
					$tally_EX_5UTR_core{$core}{$sp1}{$sp2}++;
					$combined_EX_orth_core{$core}{$sp1}{$ev}{$sp2}=1; # to get the 4-way Venn diagram 
					
                                        # just printing the positive hits
					$is_G_cl_core_5UTR{$g_cl{$sp1}{$g1}}{$sp2}{$core}=~s/\,$//;
					print "$core\t$sp1\t$sp2\tUTR5\t$ev\t$is_G_cl_core_5UTR{$g_cl{$sp1}{$g1}}{$sp2}{$core}\n";
				    }
				}
			    }
			}
		    }
		}
	    }
	}
    }
}

### 6) Prints all the stats
print "CORE\tSP1\tSP2\tTOTAL_EX\tWith_G_Orth\tWith_G_Orth_Core\tWith_Ex_Orth\tWith_Ex_Orth_Core\tEX_in_5UTR_all\tEX_in_5UTR_wG_core\tWith_5UTR_both\n";
foreach my $core (@CORES){
    foreach my $sp1 (@SPECIES){
	foreach my $sp2 (@SPECIES){
	    if ($sp1 ne $sp2){
		my $total=$#{$cores_ev{$core}{$sp1}}+1;
		$tally_G_orth_core{$core}{$sp1}{$sp2}=0 if !$tally_G_orth_core{$core}{$sp1}{$sp2};
		$tally_EX_orth{$core}{$sp1}{$sp2}=0 if !$tally_EX_orth{$core}{$sp1}{$sp2};
		$tally_EX_orth_core{$core}{$sp1}{$sp2}=0 if !$tally_EX_orth_core{$core}{$sp1}{$sp2};
		$tally_EX_5UTR_all{$core}{$sp1}=0 if !$tally_EX_5UTR_all{$core}{$sp1};
		$tally_EX_5UTR_core{$core}{$sp1}{$sp2}=0 if !$tally_EX_5UTR_core{$core}{$sp1}{$sp2};
		$tally_EX_5UTR{$core}{$sp1}=0 if !$tally_EX_5UTR{$core}{$sp1};

		print "$core\t$sp1\t$sp2\t$total\t$tally_G_orth{$core}{$sp1}{$sp2}\t$tally_G_orth_core{$core}{$sp1}{$sp2}\t".
		    "$tally_EX_orth{$core}{$sp1}{$sp2}\t$tally_EX_orth_core{$core}{$sp1}{$sp2}\t".
		    "$tally_EX_5UTR_all{$core}{$sp1}\t$tally_EX_5UTR{$core}{$sp1}\t$tally_EX_5UTR_core{$core}{$sp1}{$sp2}\n";
	    }
	}
    }
    print "\n";
}

### 7) Gets the 4-way Venn diagram
my %tally_for_gene_4way;
my %tally_for_exon_4way;
my %total_for_4way; 

foreach my $core (@CORES){
    foreach my $sp (@SPECIES){
	my @events=@{$cores_ev{$core}{$sp}};
	foreach my $ev (@events){ # One entry per exon, but the same gene could be multiple times.
	    if ($sp eq "araTha10"){
		my $spA="ce11"; my $spB="dm6"; my $spC="hg38";
		$combined_G_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringG = $combined_G_orth_core{$core}{$sp}{$ev}{$spA}.$combined_G_orth_core{$core}{$sp}{$ev}{$spB}.$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_gene_4way{$core}{$sp}{$stringG}++;
		
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringEX = $combined_EX_orth_core{$core}{$sp}{$ev}{$spA}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_exon_4way{$core}{$sp}{$stringEX}++;

		$total_for_4way{$core}{$sp}++;
	    }
	    elsif ($sp eq "ce11"){
		my $spA="araTha10"; my $spB="dm6"; my $spC="hg38";
		$combined_G_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringG = $combined_G_orth_core{$core}{$sp}{$ev}{$spA}.$combined_G_orth_core{$core}{$sp}{$ev}{$spB}.$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_gene_4way{$core}{$sp}{$stringG}++;

		$combined_EX_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringEX = $combined_EX_orth_core{$core}{$sp}{$ev}{$spA}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_exon_4way{$core}{$sp}{$stringEX}++;  

		$total_for_4way{$core}{$sp}++;
	    }
	    elsif ($sp eq "dm6"){
		my $spA="araTha10"; my $spB="ce11"; my $spC="hg38";
		$combined_G_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringG = $combined_G_orth_core{$core}{$sp}{$ev}{$spA}.$combined_G_orth_core{$core}{$sp}{$ev}{$spB}.$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_gene_4way{$core}{$sp}{$stringG}++;

		$combined_EX_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringEX = $combined_EX_orth_core{$core}{$sp}{$ev}{$spA}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_exon_4way{$core}{$sp}{$stringEX}++;

		$total_for_4way{$core}{$sp}++;
	    }
	    elsif ($sp eq "hg38"){
		my $spA="araTha10"; my $spB="ce11"; my $spC="dm6";
		$combined_G_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_G_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringG = $combined_G_orth_core{$core}{$sp}{$ev}{$spA}.$combined_G_orth_core{$core}{$sp}{$ev}{$spB}.$combined_G_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_gene_4way{$core}{$sp}{$stringG}++;

		$combined_EX_orth_core{$core}{$sp}{$ev}{$spA}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spA};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spB};
		$combined_EX_orth_core{$core}{$sp}{$ev}{$spC}=0 if !$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		my $stringEX = $combined_EX_orth_core{$core}{$sp}{$ev}{$spA}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spB}.$combined_EX_orth_core{$core}{$sp}{$ev}{$spC};
		$tally_for_exon_4way{$core}{$sp}{$stringEX}++;

		$total_for_4way{$core}{$sp}++;
	    }
	}
    }
}

my @STRINGS = ("000","100","010","001","110","101","011","111");
print "CORE\tSPECIES\tCONS_STRING\tTOTAL_EV\tW_G_Orth\tPerc_W_G_Orth\tW_EX_Orth\tPerc_W_EX_Orth\n";
foreach my $core (@CORES){
    foreach my $sp (@SPECIES){
	foreach my $string (@STRINGS){
	    $tally_for_gene_4way{$core}{$sp}{$string}=0 if !$tally_for_gene_4way{$core}{$sp}{$string};
	    $tally_for_exon_4way{$core}{$sp}{$string}=0 if !$tally_for_exon_4way{$core}{$sp}{$string};
	    my $perc_G = sprintf("%.2f",100*$tally_for_gene_4way{$core}{$sp}{$string}/$total_for_4way{$core}{$sp});
	    my $perc_EX = sprintf("%.2f",100*$tally_for_exon_4way{$core}{$sp}{$string}/$total_for_4way{$core}{$sp});
	    print "$core\t$sp\t$string\t$total_for_4way{$core}{$sp}\t$tally_for_gene_4way{$core}{$sp}{$string}\t$perc_G\t$tally_for_exon_4way{$core}{$sp}{$string}\t$perc_EX\n";
	}
    }
    print "\n";
}
