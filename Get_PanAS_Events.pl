#!/usr/bin/perl
# Original name: DistribAS-BINS-byPSI_v3.pl
use strict;
no strict "refs";
use warnings;
use Getopt::Long;

### General options & defaults
my $sp="Ath";
my $min=10;
my $max=90;
my $min_samples=20;
my $Q="O[WK]\,.+?\,.+?\,.+?\,.+?\@"; # VTS coverage score
my $typeFILE="AltEx";
my $helpFlag;
my $no_pIR;
my $VLOW;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "event_type=s" => \$typeFILE,
			  "min_samples=i" => \$min_samples,
                          "help" => \$helpFlag,
                          "no_pIR" => \$no_pIR,
                          "VLOW" => \$VLOW
    );

if (!defined($ARGV[0]) || $helpFlag){
    die "\nUsage: DistribAS-BINS.pl INCLUSION [options]

Identifies PanAS, MildAS and Switch AS events as per Tapial et al 2017

[General options] 
        --event_type s           Type of AS event to be tested: AltEx, IR, Alt3, Alt5 (default AltEx)
        --VLOW                   Allows using samples with VLOW coverage (default OFF)
        --no_pIR                 Does NOT apply the binomial test filter for IR (default OFF)

*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

### Just to get the GeneID form the VastDB file
open (KEY, "Ath.Event-Gene.IDs.txt") || die "Needs the Event-Gene key\n";
my %key;
while (<KEY>){
    chomp;
    my @t=split(/\t/);
    $key{$t[0]}=$t[1];
}
close KEY;

open (I, $ARGV[0]) || die "Needs an INCLUSION_LEVELS Table\n";
my $head=<I>;
chomp($head);
my @head=split(/\t/,$head);
my $tis=join("\t",@head[6..$#head]);
my %Q; my %val; # conversion from Sample name to array index (too complex for historical reasons)
for my $i (6..$#head){
    if ($head[$i]=~/(.+?)\-Q/){
        $Q{$1}=$i;
    }
    else {
        if ($head[$i]=~/(.+?)\-C/){$head[$i]=$1;}
	$val{$head[$i]}=$i;
    }
}

### Output files
open (DATA, ">Events_$sp-$min_samples-$typeFILE-DATA-PastDB.tab");
open (ALL_DATA, ">Events_$sp-$min_samples-$typeFILE-ALL_DATA-PastDB.tab");

open (EV, ">EVENTS_PASTDB/Events_$sp-$min_samples-80-$typeFILE-PastDB.txt");
open (EV2, ">EVENTS_PASTDB/Events_$sp-$min_samples-10_25-$typeFILE-PastDB.txt");
open (EVswitch, ">EVENTS_PASTDB/Events_$sp-$min_samples-Switch-$typeFILE-PastDB.txt");
open (EVBG, ">EVENTS_PASTDB/Events_$sp-$min_samples-BG-$typeFILE-PastDB.txt");

open (G, ">GENES_PASTDB/Genes_$sp-$min_samples-80-$typeFILE-PastDB.txt");
open (G2, ">GENES_PASTDB/Genes_$sp-$min_samples-10_25-$typeFILE-PastDB.txt");
open (Gswitch, ">GENES_PASTDB/Genes_$sp-$min_samples-Switch-$typeFILE-PastDB.txt");
open (GBG, ">GENES_PASTDB/Genes_$sp-$min_samples-BG-$typeFILE-PastDB.txt") || die "BG genes";

print DATA "GENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tN_samples\tPercAS\tAvPSI\tMIN\tMAX\tRANGE\tTYPE\n";
print ALL_DATA "GENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tN_samples\tPercAS\tAvPSI\tMIN\tMAX\tRANGE\n";
#######

### Starts processing the INCLUSION table
my %name; my %ev_co; my %fullco; # Event basic info storage
my %PSI_data; # stores PSIs for the event
my (%sumT, %data);
my %tally; # count of the of % AS samples
my (%doneG, %doneG2, %doneGsw, %doneGBG); # to avoid printing out genes multiple times

while (<I>){
    s/VLOW/N/g unless $VLOW;
    chomp;
    my @t=split(/\t/);
    my $type;

    if ($t[1]=~/EX/){
	$type="AltEx";
    }
    elsif ($t[1]=~/INT/){
	$type="IR";
    }
    else {
	$type=$t[5];
    }
    
    my $ev=$t[1];
    $name{$ev}=$t[0];
    $ev_co{$ev}=$t[2];
    $fullco{$ev}=$t[4];
    my $g;
    $g=$key{$ev};
    next if !$g;

    my ($maxEV,$minEV); # max and min PSI values the event takes
    my ($sumEV, $AS, $noAS, $tot, $ASmaxi) = (0,0,0,0,0); # different counters for the event
    my $pre=join("\t",@t[0..5]);
    
    foreach my $tis (sort (keys %val)){ # complex for historial reasons
	if ($t[$Q{$tis}]=~/$Q/){
	    my $kill_pIR = "";
	    if ($t[5]=~/IR/ && !$no_pIR){
                my ($temp_p)=$t[$Q{$tis}]=~/\,.+?\,.+?\,.+?\,(.+?)\@/;
                $kill_pIR = 1 if $temp_p < 0.05;
            }
	    next if $kill_pIR;

	    $tot++; # number of samples with coverage

	    if ($type =~ /$typeFILE/){
		$AS++ if ($t[$val{$tis}]>=$min && $t[$val{$tis}]<=$max); #not counting switches
		$noAS++ if ($t[$val{$tis}]<=$min || $t[$val{$tis}]>=$max);
		
		$sumEV+=$t[$val{$tis}];		
                if (defined $maxEV){
                     $maxEV=$t[$val{$tis}] if $t[$val{$tis}]>=$maxEV;
                } else {
                     $maxEV=$t[$val{$tis}];
                }
                if (defined $minEV){
		     $minEV=$t[$val{$tis}] if $t[$val{$tis}]<=$minEV;
                } else {
                     $minEV=$t[$val{$tis}];
                }
		$PSI_data{$t[1]}.="$t[$val{$tis}]\n";
	    }
	}
    }
    if ($tot >= $min_samples){
	if ($type =~ /$typeFILE/){
	    my $perc = sprintf("%.0f",100*$AS/$tot);
	    $tally{$perc}++;
	    my $avEV = sprintf("%.2f",$sumEV/$tot);
	    $sumT{$perc}+=$avEV;
	    $data{$perc}.=$PSI_data{$t[1]};
	    
	    # for Switch
	    my $perc_noAS=sprintf("%.0f",100*$noAS/$tot);
	    my $range=$maxEV-$minEV;
	    
	    if ($perc > 80){ # > to match the bins
		print DATA "$pre\t$tot\t$perc\t$avEV\t$minEV\t$maxEV\t$range\tPAN_AS_80\n";
		print EV "$ev\n";
		print G "$g\n" if !$doneG{$g};
		$doneG{$g}=1;
	    }
	    elsif ($perc > 10 && $perc <=25){ # 10% of the samples to adjust it to MIC paper AS def
		print DATA "$pre\t$tot\t$perc\t$avEV\t$minEV\t$maxEV\t$range\t10_25\n";
		print EV2 "$ev\n";
		print G2 "$g\n" if !$doneG2{$g};
		$doneG2{$g}=1;
	    }
	    # switch 
	    if ($perc_noAS > 80 && $range > 80) {
		print DATA "$pre\t$tot\t$perc\t$avEV\t$minEV\t$maxEV\t$range\tSWITCH\n";
		print EVswitch "$ev\n";
		print Gswitch "$g\n" if !$doneGsw{$g};
		$doneGsw{$g}=1;
	    }
	    
	    print ALL_DATA "$pre\t$tot\t$perc\t$avEV\t$minEV\t$maxEV\t$range\n";
	    print EVBG "$ev\n"; # any event with matched coverge filter
	}
	print GBG "$g\n" if !$doneGBG{$g};
	$doneGBG{$g}=1;
    }
}

my ($option_VLOW,$option_pIR);
$option_VLOW="YES" if $VLOW;
$option_VLOW="NO" if !$VLOW;
$option_pIR="YES" if !$no_pIR;
$option_pIR="NO" if $no_pIR;

print "$typeFILE; >=$min_samples samples; $min-$max; VLOW: $option_VLOW; p_IR: $option_pIR\n";
print "FRACTION\tN\tAVERAGE\n";
my $sum=0;
my $sumAV=0;
my $dataEV;
for my $i (1..100){
    if ($i%5==0){
	$sum+=$tally{$i};
	$sumAV+=$sumT{$i};
	$dataEV.=$data{$i};
	my $p=$i-4;
	my $head="$p-$i";
	my $AV=sprintf("%.2f",$sumAV/$sum);
	print "$head\t$sum\t$AV\n";

	open ($head,">DATA_PASTDB/$sp-$typeFILE-$min_samples-DATA_PSI-$head.txt");
	print $head "$dataEV";
	close $head;

	$dataEV="";
	$sumAV=$sum=0;
    }
    else {
	$sum+=$tally{$i};
	$sumAV+=$sumT{$i};
	$dataEV.=$data{$i};
    }
}
