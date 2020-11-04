#!/usr/bin/perl
# Original name: TrulyAS2.pl
use strict;
use warnings;
use Getopt::Long;

### Setting default variables
my $range_PSI = 25; #for trulyAS
my $min_PSI = 10;
my $max_PSI = 90;
my $fract_AS = 0.1;
my $N_expr = 3;
my $VLOW;
my $no_pIR;
my $helpFlag;
my $min_ALT_use = 0;
my $command = join(" ", @ARGV);
my $outFile;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(  "range_PSI=i"   => \$range_PSI,
	     "min_PSI=i"     => \$min_PSI,
	     "max_PSI=i"     => \$max_PSI,
	     "fract_AS=f"    => \$fract_AS,
	     "N=i"           => \$N_expr,
	     "VLOW"          => \$VLOW,
	     "no_pIR"        => \$no_pIR,
	     "min_ALT_use=i" => \$min_ALT_use,
	     "outFile=s"     => \$outFile,
	     "help"          => \$helpFlag,
	     "h"             => \$helpFlag
    );


die "
Usage: Get_Event_Stats.pl INCLUSION_TABLE.tab [options]

OPTIONS:

       --N i                Minimum number of samples with coverage (def = 3)
       --VLOW               Use VLOW coverage (def = OFF)
       --no_pIR             Do not filter IR by binomial p-value (def = OFF)
       --range_PSI i        Minimum range to be considered TrulyAS (SuperAS is 2x) (def = 25)
       --fract_AS fr        Fraction of samples with AS to be TrulyAS (def = 0.1) 
       --min_PSI i          Low threshold for an event to be TrulyAS (SuperAS is 2x) (def = 10)
       --max_PSI i          High threshold for an event to be TrulyAS 
                                 (SuperAS is 100-2*[100-max_PSI]) (def = 90)
       --min_ALT_use i      Minimum inclusion of the exon in which the Alt3/Alt5 is located across all 
                                 compared samples (default 0) (combine >= v2.2.1)
       --outFile            File name (if not provided, constructed by default)
       -h, --help           Print this help message


* Note: Default options for VLOW and p_IR options are the opposite as in vast-tools

* Definitions of scores (in all cases, for events with >= N valid samples):
    - SuperAS = if fract_AS of samples have a min_PSI*2 <= PSI <= 100-2*[100-max_PSI] 
                   OR a range of PSIs >= 2*range_PSI. 
    - TrulyAS = if fract_AS of samples have a min_PSI <= PSI <= max_PSI 
                   OR a range of PSIs >= range_PSI and are not SuperAS.
    - HIGH_PSI = if min_PSI > 90.
    - LOW_PSI = if max_PSI < 10.
    - MILD_AS_H = else, if average PSI > 50.
    - MILD_AS_L = else, if average PSI < 50.

*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

" if !$ARGV[0] || $helpFlag;

open (I, $ARGV[0]) || die "Cannot open input file as ARGV[0]\n";

unless ($outFile){
    my ($root)=$ARGV[0]=~/([^\/]+)\.t[ax][bt]/;
    $root.="-noVLOW" if !$VLOW;
    $root.="-pIR" if !$no_pIR;
    $root.="-min_ALT_use$min_ALT_use" if $min_ALT_use > 0;
    $root.="-N$N_expr";
    
    $outFile = "$root.trulyAS2";
}
open (O, ">$outFile") || die "Cannot open output file $outFile\n";

my $head=<I>;
chomp($head);
my @head=split(/\t/,$head);
my %Q; my %val;
for my $i (6..$#head){
    if ($head[$i]=~/(.+?)\-Q/){
	$Q{$1}=$i;
    }
    else {
	if ($head[$i]=~/(.+?)\-C/){$head[$i]=$1;}
	$val{$head[$i]}=$i;
    }
}

my %tally;
my %all_types;
my %tally_str;

print O "GENE\tEVENT\tCOORD\tLENGTH\tFullCO\tCOMPLEX\tTotSamples\tN_AS\tN_SAS\tTYPE\tAverage\tMin\tMax\tRange\n";
while (<I>){
    s/VLOW/N/g if !$VLOW;
    chomp;
    my @t=split(/\t/);
    my $ev=$t[1];
    my $pre=join("\t",@t[0..5]); # basic info
    
    my ($SAS, $sum, $AS, $tot) = (0,0,0,0);
    my ($min,$max);
    foreach my $tis (sort (keys %val)){
	### Valid for old and new coverage scores
	if (($t[$Q{$tis}]=~/O[KW]\,.+?\,.+?\,.+?\,/ && $t[$Q{$tis}]!~/\@/) || ($t[$Q{$tis}]=~/O[KW]\,.+?\,.+?\,.+?\,.+?\@/ && $t[$Q{$tis}]=~/\@/)){
	    my $kill_pIR = "";
	    if ($t[5]=~/IR/ && !$no_pIR){
		my ($temp_p)=$t[$Q{$tis}]=~/O[KW]\,.+?\,.+?\,.+?\,(.+?)\@/;
		$kill_pIR = 1 if $temp_p < 0.05;
	    }

	    my $kill_ALT = "";
	    if ($t[5]=~/Alt[35]/ && $min_ALT_use > 0){
		my ($temp_ALT)=$t[$Q{$tis}]=~/O[KW]\,.+?\,(.+?)\,.+?\,.+?\@/;
		if ($temp_ALT=~/\d/){
		    $kill_ALT = 1 if $temp_ALT < $min_ALT_use;
		} 
		# else => old INCL version
	    }
	    
	    unless ($kill_pIR || $kill_ALT){
		$tot++;
		$AS++ if ($t[$val{$tis}]>=$min_PSI && $t[$val{$tis}]<=$max_PSI);
		$SAS++ if ($t[$val{$tis}]>=$min_PSI+$min_PSI && $t[$val{$tis}]<=$max_PSI-$min_PSI);
		
		if (defined $min){
		    $min=sprintf("%.2f",$t[$val{$tis}]) if $t[$val{$tis}]<=$min;
		} else {
		    $min=sprintf("%.2f",$t[$val{$tis}]);
		}
		if (defined $max){
		    $max=sprintf("%.2f",$t[$val{$tis}]) if $t[$val{$tis}]>=$max;
		} else {
		    $max=sprintf("%.2f",$t[$val{$tis}]);
		}
		$sum+=$t[$val{$tis}];
	    }
	}
    }

    if ($tot >= $N_expr){
	my $range = sprintf("%.2f",$max-$min);
	my $av = sprintf("%.2f",$sum/$tot);
	my $type=""; my $type_ev=""; my $str="";

	# gets strand for quality purposes
	if ($t[5]=~/IR/){
	    ($str)=$t[4]=~/.+\:([\+\-])$/;
	    $type_ev="IR";
	}
	elsif ($t[5]=~/Alt3/){
	    #chr3R:16827135,16826051-16826307+16826304
	    my ($C1d)=$t[4]=~/\:(\d+)/;
	    my $C2a;
	    ($C2a)=$t[4]=~/\,(\d+)/;
	    ($C2a)=$t[4]=~/\-(\d+)/ if !$C2a;
	    $str="+" if $C1d<$C2a;
	    $str="-" if $C1d>$C2a;
	    $type_ev="Alt3";
	}
	elsif ($t[5]=~/Alt5/){
	    #chr2R:22172082-22172248+22172252,22172316
	    my ($C1d,$C2a);
	    ($C1d)=$t[4]=~/\:(\d+)/;
	    ($C1d)=$t[4]=~/\-(\d+)/ if !$C1d;
	    ($C2a)=$t[4]=~/\,(\d+)/;
	    $str="+" if $C1d<$C2a;
	    $str="-" if $C1d>$C2a;
	    $type_ev="Alt5";
	}
	else {
	    my ($C1d,$C2a)=$t[4]=~/(\d+)\,.+\,(\d+)/;
	    $str="+" if $C1d<$C2a;
	    $str="-" if $C1d>$C2a;
	    $type_ev="AltEx";
	}

	if ($SAS/$tot>=$fract_AS || $range>=$range_PSI*2){$type="SuperAS"; $tally_str{$str}++;} #>10% are SuperAS or range>50
	elsif ($AS/$tot>=$fract_AS || $range>=$range_PSI){$type="TrulyAS";}
	elsif ($max<10){$type="LOW_PSI";}
	elsif ($min>90){$type="HIGH_PSI";}
	else {
	    $type="MILD_AS_L" if $av<50;
	    $type="MILD_AS_H" if $av>=50;
	}
	
	$tally{$type_ev}{$type}++; # modified on 07/03/16
	$all_types{$type}=1;

	print O "$pre\t$tot\t$AS\t$SAS\t$type\t$av\t$min\t$max\t$range\n";
    }
}

my ($using_VLOW, $using_pIR);
$using_VLOW = "NO";
$using_pIR = "NO";
$using_VLOW = "YES" if $VLOW;
$using_pIR = "YES" if !$no_pIR;

print "Command: TrulyAS2.pl $command

Options: 
- Min samples coverage (N):   $N_expr
- Min fraction for Truly:     $fract_AS
- Min range of PSIs:          $range_PSI
- Low bound PSI for AS:       $min_PSI
- High bound PSI for AS:      $max_PSI
- Using VLOW scores?          $using_VLOW
- Filtering IR by p_IR?       $using_pIR
- Min ALT PSI required        $min_ALT_use

";

print ">>> Events by type:\n";
foreach my $type_ev (sort {$b cmp $a} (keys %tally)){
    print "  - $type_ev:\n";
    foreach my $type (sort {$b cmp $a} (keys %all_types)){
	$tally{$type_ev}{$type}=0 if !$tally{$type_ev}{$type};
	print "    $type\t$tally{$type_ev}{$type}\n";
    }
    print "\n";
}

print "\n>>> SAS strand-check: (+) $tally_str{'+'}  (-) $tally_str{'-'}\n\n";

