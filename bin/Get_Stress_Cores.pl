#!/usr/bin/perl
use strict;
use warnings;

### Input file: a processed vast-tools table for the 5 experiments:
# Columns:
# 1: GeneID
# 2-6: standard 6 columns from vast-tools INCLUSION TABLE
# 7-11: dPSIs for each experiment (5). NA = no coverage in control and/or stress
# 12-21: PSIs for each sample (10). NA = no coverage.
my $input = $ARGV[0];
my ($stress_type) = $input =~ /(.+?)\-dPSI_PSI-table.tab/;

die "
Usage: Get_AS_BG_etc_Abiotic.pl STRESS-dPSI-PSI_Table.tab

It uses a processed table as input: STRESS_dPSI-PSI_Table  
 # Columns:
 # 1: GeneID
 # 2-6: standard 6 columns from vast-tools INCLUSION TABLE
 # 7-11: dPSIs for each experiment (5). NA = no coverage in control and/or stress
 # 12-21: PSIs for each sample (10). NA = no coverage.

Cut-offs defined inside the script.

" if !$ARGV[0];

# General cut-offs
my $max_dPSI=5; # to be considered non-change
my $min_PSI=10; my $max_PSI=90; # to be considered as "AS event"
my $min_cov=2; # min number of experiments for which there must be coverage
my $min_AS_exp=1; # min number of experiments in which the event has to be AS

# Cut-offs to define the core set:
my $core_cutoff = 15; # dPSI to be considered regulated
my $core_minN = 2; # min number of experiments in which is regulated

# Prepares the output file name
my $output_file = "$stress_type"."_core-AS_NC-min_cov$min_cov-min_dPSI$max_dPSI-AS$min_PSI"."_$max_PSI-min_AS_exp$min_AS_exp-Core_minN$core_minN-Core_cutoff$core_cutoff.tab";
open (O, ">$output_file") || die "It cannot open output";
print O "GENE_ID\tGENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tN_cov\tChange_UP\tChange_DOWN\tN_low_dPSI\tN_AS_exp\tN_AS_samples\tTYPE\n";

my %tally; # tally of event by type
open (I, $input) || die "Needs file ($input)\n";
<I>;
while (<I>){
    s/\r//g; # if the input was generated using excel
    s/\"//g;
    chomp;
    my @t=split(/\t/);

    my $event=$t[2];
    next if !$event; # skips empty lines

    my ($type)=$event=~/Ath([^\d]+?)\d/;
    
    ### Initializes variables of interest per event
    my $N_cov = 0; # number of exp with coverage
    my $N_low_dPSI = 0; # number of expr with dPSI < $max_dPSI
    my $N_AS_samples = 0; # number of samples in which it is AS
    my $N_AS_exp = 0; # number of exp in which it is AS
    my ($core_UP,$core_DOWN)=(0,0); # count of differentially spliced UP or DOWN

    ### first, checks dPSI and coverage
    for my $i (7..11){ # loops through the 5 dPSI values
	if ($t[$i]=~/\d/){ # i.e. checks if it has a value (number = \d) (i.e. non-NA)
	    $N_cov++;
	    
	    if (abs($t[$i])<$max_dPSI){ # checks if absolute dPSI < 10
		$N_low_dPSI++;
	    }
	    elsif ($t[$i] >= $core_cutoff){
		$core_UP++;
	    }
	    elsif ($t[$i] <= -$core_cutoff){
		$core_DOWN++;
	    }
	}
    }
    
    ### Now defines which are AS based on PSI levels
    for my $i (12..21){
	if ($t[$i]=~/\d/){
	    $N_AS_samples++ if ($t[$i]>$min_PSI && $t[$i]<$max_PSI); # i.e. it's between 10 and 90
	}
	# Now it analyzes it by pairs
	if ($i%2==0){ # i.e. $i is multiple of 2, i.e. the control
	    ### Adds one if $i or $i+1 (cont or stress) are AS
	    if ($t[$i]=~/\d/ && $t[$i+1]=~/\d/){
		$N_AS_exp++ if (($t[$i]>$min_PSI && $t[$i]<$max_PSI) || ($t[$i+1]>$min_PSI && $t[$i+1]<$max_PSI));
	    }
	}
    }
    
    ### puts stuff together
    my $regulation="NO_COV";
    my $line=join("\t",@t[0..6]);
    
    if ($N_cov>= $min_cov){
	$regulation="OTHER_COV";
	### is AS_noCH if: it has at least $min_cov (2) experiments with coverage,
	###                in all of them the abs dPSI is < $max_dPSI (10)
	###                it's AS in at least $min_AS_exp (2) experiments 
	$regulation="AS_NO_CH" if $N_cov>= $min_cov && $N_AS_exp >= $min_AS_exp && $N_low_dPSI == $N_cov;

	$regulation="CORE_UP" if $core_UP >= $core_minN && $core_DOWN == 0;
	$regulation="CORE_DOWN" if $core_DOWN >= $core_minN && $core_UP == 0;
	$regulation="AMBIGUOUS" if ($core_UP>0 && $core_DOWN>0);
    }
    print O "$line\t$N_cov\t$core_UP\t$core_DOWN\t$N_low_dPSI\t$N_AS_exp\t$N_AS_samples\t$regulation\n";
    $tally{$regulation}{$type}++;
}

my @TYPES=("EX","INT","ALTD","ALTA");
my @REG=("CORE_UP","CORE_DOWN","AMBIGUOUS","AS_NO_CH","OTHER_COV","NO_COV");

print "SETTINGS\n".
    "\tCOVERAGE: Min_cov = $min_cov\n".
    "\tCORE: Min_N_changing = $core_minN, dPSI cutoff = 15\n".
    "\tAS_NO_CH: Min_AS_exp = $min_AS_exp, Max_dPSI = $max_dPSI, PSI = $min_PSI-$max_PSI\n\n";
print "RESULTS ($stress_type core)\n";
print "TYPE\tEX\tINT\tALTD\tALTA\tTOTAL\n";
foreach my $reg (@REG){
    my $total=0;
    print "$reg";
    
    foreach my $type (@TYPES){
	$total+=$tally{$reg}{$type};
	print "\t$tally{$reg}{$type}";
    }
    print "\t$total\n";
}
