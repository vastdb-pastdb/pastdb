#!/usr/bin/perl
# Original name: QuantifyAS_byPSI_v2.pl
use strict;
use warnings;
use Getopt::Long;

my $event_gene=$ARGV[1]; # Ath.Event-Gene.IDs.txt

### Default variables
my $test_type_AS = "ALL";
my $step=1;
my $rep=100;
my $noVLOW;
my $no_pIR;
my $helpFlag;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(   "AS_type=s" => \$test_type_AS,
	      "step_size=i"    => \$step,
	      "iter=i"         => \$rep,
	      "noVLOW"         => \$noVLOW,
	      "no_pIR"         => \$no_pIR,
	      "help"           => \$helpFlag,
	      "h"              => \$helpFlag
    );


if (!$event_gene || defined ($helpFlag)){
    die "
Usage: Quantify_AS_by_Subsampling.pl INCLUSION_FULL_TABLE Event_Gene [Options]

OPTIONS:
       --AS_type XX            Type of AS to be tested: ALL, EX, INT, ALTA, ALTD (def = ALL)
       --step_size i           Step increase in number of samples (def = 1) 
       --iter i                Number of iterations (def = 100)
       --noVLOW                Excludes VLOW coverage scores 
       --no_pIR                It does NOT perform the imbalance test for IR

";
}

open (I, $ARGV[0]) || die "INCLUSION TABLE\n";
my ($root)=$ARGV[0]=~/.+\-(.+?)[\-\_\.]/;

open (G, $event_gene) || die "Needs the Event to Gene ID\n";
my %ev_g;
while (<G>){
    chomp;
    my @t=split(/\t/);
    $ev_g{$t[0]}=$t[1];
}
close G;


my %PSI; my %tot; my @valid_ev;
my $last_sample_number;
<I>;
while (<I>){
    chomp;
    s/VLOW/N/g if $noVLOW;
    my @t=split(/\t/);
    my $ev=$t[1];

    for my $i (6..$#t){
	if ($i%2 == 0){
	    my $sample=($i-6)/2;
	    $last_sample_number = $sample;
	    if ($t[$i+1]=~/O[KW]\,.+?\,.+?\,.+?\,.+?\@/){
		my $kill_pIR = "";
		if ($t[5]=~/IR/ && !$no_pIR){
		    my ($temp_p)=$t[$i+1]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
		    $kill_pIR = 1 if $temp_p < 0.05;
		}
		
		unless ($kill_pIR){
		    $PSI{$ev}[$sample]=$t[$i];
		    $tot{$ev}++;
		}
		else {
		    $PSI{$ev}[$sample]="NA";
		}
	    }
	    else {
		$PSI{$ev}[$sample]="NA";
	    }
	}
    }
    if (defined $tot{$ev}){
	push(@valid_ev,$ev) if $tot{$ev}>=3;
    }
}
close I;

### gets the master array with samples
my @master;
for my $a (0..$last_sample_number){
    $master[$a]=$a;
}

my $n_samples=$last_sample_number;
my $limit=30; # artificial upper bound (30)

if ($noVLOW){
    open (O1, ">Perc_AS_genes_bins-$root-noVLOW-Step_$step-Iter_$rep-MinAS_10-90-$test_type_AS.txt");
    open (O2, ">Perc_AS_genes_bins-$root-noVLOW-Step_$step-Iter_$rep-MinAS_20-80-$test_type_AS.txt");
    open (O3, ">Perc_AS_genes_bins-$root-noVLOW-Step_$step-Iter_$rep-MinAS_30-70-$test_type_AS.txt");
}
else {
    open (O1, ">Perc_AS_genes_bins-$root-Step_$step-Iter_$rep-MinAS_10-90-$test_type_AS.txt");
    open (O2, ">Perc_AS_genes_bins-$root-Step_$step-Iter_$rep-MinAS_20-80-$test_type_AS.txt");
    open (O3, ">Perc_AS_genes_bins-$root-Step_$step-Iter_$rep-MinAS_30-70-$test_type_AS.txt");
}
print O1 "N_samples\tPerc_AS_genes\n";
print O2 "N_samples\tPerc_AS_genes\n";
print O3 "N_samples\tPerc_AS_genes\n";

for (my $i = 5; $i <= $limit; $i+=$step){
    for my $j (1..$rep){
	print "Doing Step $i Iteration $j\n";
	my @temp=@master;
	my @selected=();
	# now it gets $i random samples
	for my $k (1..$i){
	    my $index=int(rand($#temp+1));
	    push(@selected,$temp[$index]);
	    splice(@temp,$index,1);
	}
	
	my %doneG=(); my $count_valid_g=0;
	my %doneAS1=(); my $count_AS_g1=0;
	my %doneAS2=(); my $count_AS_g2=0;
	my %doneAS3=(); my $count_AS_g3=0;

	foreach my $ev (@valid_ev){
	    my ($type_AS)=$ev=~/^.{3}(.+?)\d/;
	    # for each event
	    my $tot_cov=0;
	    my ($tot_AS1, $tot_AS2, $tot_AS3)=(0,0,0);
	    
	    my $g=$ev_g{$ev};
	    next if !$g;
	    
	    foreach my $sample (@selected){
		if ($PSI{$ev}[$sample] ne "NA"){
		    $tot_cov++;
		    
		    if ($test_type_AS eq "ALL" || $type_AS eq $test_type_AS){
			$tot_AS1++ if ($PSI{$ev}[$sample]>=10 && $PSI{$ev}[$sample]<=90);
			$tot_AS2++ if ($PSI{$ev}[$sample]>=20 && $PSI{$ev}[$sample]<=80);
			$tot_AS3++ if ($PSI{$ev}[$sample]>=30 && $PSI{$ev}[$sample]<=70);
		    }
		}
	    }
	    if ($tot_cov >= 5){
		### AS Gene
		if ($tot_AS1 > 0){
		    $count_AS_g1++ if !$doneAS1{$g};
		    $doneAS1{$g}=1;
		}
		if ($tot_AS2 > 0){
		    $count_AS_g2++ if !$doneAS2{$g};
		    $doneAS2{$g}=1;
		}
		if ($tot_AS3 > 0){
		    $count_AS_g3++ if !$doneAS3{$g};
		    $doneAS3{$g}=1;
		}

		### Valid Gene
		$count_valid_g++ unless $doneG{$g};
		$doneG{$g}=1;
	    }
	}
	
	### Now gets the % of AS genes
	my $perc_AS_genes1 = sprintf("%.2f",100*$count_AS_g1/$count_valid_g);
	my $perc_AS_genes2 = sprintf("%.2f",100*$count_AS_g2/$count_valid_g);
	my $perc_AS_genes3 = sprintf("%.2f",100*$count_AS_g3/$count_valid_g);

	print O1 "s$i\t$perc_AS_genes1\n";
	print O2 "s$i\t$perc_AS_genes2\n";
	print O3 "s$i\t$perc_AS_genes3\n";
    }
}
