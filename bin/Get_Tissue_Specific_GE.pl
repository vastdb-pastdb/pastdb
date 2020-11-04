#!/usr/bin/perl
# original name: get_tisGE_v3.pl
# group file: Sample_name\tGroup\tGr_Excluded1,Gr_Excluded2

use strict;
no strict "refs";
use warnings;
use Getopt::Long;

#### General variables
my $min_expr=5; # to define minimal GE of maximum
my $min_fold_glob=5;
my $min_fold=3;
my $min_range=2; # specially for DOWN
my $min_rep=1;
my $print;
my $average;
my $test_tis;
my $groups;
my $strict;
my $helpFlag;

Getopt::Long::Configure("no_auto_abbrev");
GetOptions(               "min_fold=i" => \$min_fold,
			  "min_fold_glob=i" => \$min_fold_glob,
                          "min_expr=i" => \$min_expr,
			  "min_range=i" => \$min_range,
                          "groups=s" => \$groups,
                          "g=s" => \$groups,
			  "print" => \$print,
			  "strict" => \$strict,
			  "average" => \$average,
			  "test_tis=s" => \$test_tis,
			  "min_rep=i" => \$min_rep,
                          "help" => \$helpFlag
    );

if (!defined($ARGV[0]) || $helpFlag){
    die "\nUsage: Get_Tissue_Specific_GE.pl cRPKM-SpV-NORM.tab -g Tissue_groups [options]

Identifies tissue-specific genes by assessing non-overlapping expression among groups

[General options] 
        --min_expr i             Minimum cRPKM for the gene to be considered (default 5)
        --min_fold i             Minimum fold change of the median GE against each of the other medians (default 3)
        --min_fold_glob i        Minimum fold change of the median GE against the median of the other medians (default 5)
        --min_range i            Minimum absolute cRPKM difference with every other tissue (default 2)
        --groups/-g file         Groups file (if not provided, each sample is a group)
        --test_tis               Comma-separated list of issues for which it will make a table with fold changes (default OFF)
        --strict                 Requires all sub-samples to also be non-overlapping within fold GE 2 (default OFF)
        --average                Uses cRPKM average of group, instead of median (default OFF)
        --min_rep                Minimum number of replicates with coverage for a group to be considered (default 1)
        --print                  Print an output file with the tissue-specific genes (default OFF)


*** Questions \& Bug Reports: Manuel Irimia (mirimia\@gmail.com)

";
}

my ($v)=$ARGV[0]=~/.+?\-(.{3}.+?)[\-\.]/; # to allow subsets of V
my ($root)=$ARGV[0]=~/(.+)\./;
my ($sp)=$v=~/^(.{3})/;

print "*** Version: $v; Root: $root\n";

$root=~s/.+\///;
my $output_file="TS_GE-$v-FC$min_fold-Min_Expr$min_expr-Min_range$min_range-Rep$min_rep";
$output_file.="-Average" if $average;
$output_file.="-Strict" if $strict;

if (defined $print){
    open (O, ">$output_file-PastDB.tab");
    if (defined $average){
	print O "GENE_ID\tGENE_NAME\tTissue\tDirection\tAv_Tis\tAv_Others\tlog2FC_Tis-Others\tSD_Tis\tSD_Others\tMin_GE_Others\tMax_GE_Others\n";
    } else {
	print O "GENE_ID\tGENE_NAME\tTissue\tDirection\tMed_Tis\tMed_Others\tlog2FC_Tis-Others\tSD_Tis\tSD_Others\tMin_GE_Others\tMax_GE_Others\n";
    }
}

### Creates the files for each test tissues
my @TEST_TIS=split(/\,/,$test_tis);
my %test_tis_OK;
if ($test_tis){
    foreach my $tis (@TEST_TIS){
	my $handler="FC_$tis";
	if (defined $average){
	    open ($handler, ">FC_$tis-$v-Rep$min_rep-Average-PastDB.tab");
	    print $handler "GENE_ID\tGENE_NAME\tAv_$tis\tAv_Others\tlog2FC_$tis-Others\tSD_"."$tis\tSD_Others\tMin_Expr_Others\tMax_Expr_Others\n";
	} else {
	    open ($handler, ">FC_$tis-$v-Rep$min_rep-Median-PastDB.tab");
	    print $handler "GENE_ID\tGENE_NAME\tMed_$tis\tMed_Others\tlog2FC_$tis-Others\tSD_"."$tis\tSD_Others\tMin_Expr_Others\tMax_Expr_Others\n" if !$average;
	}    
	$test_tis_OK{$tis}=1;
    }
}

my %group;
my %temp_excl; my %excl;
if ($groups){
    open (GROUPS, $groups) || die "Can't open group file";
    while (<GROUPS>){
	s/\r//g;
	s/\"//g;
	s/ //g;
	chomp;
	my @t=split(/\t/);
	$group{$t[0]}=$t[1] if $t[1];
	$group{$t[0]}=$t[0] if !$t[1];
	
	if (defined $t[2]){
            if (defined $group{$t[0]}){
		if (defined $temp_excl{$group{$t[0]}}){
		    die "Inconsistent exclusion set for the same group ($group{$t[0]}: $temp_excl{$group{$t[0]}})\n" if $temp_excl{$group{$t[0]}} ne $t[2] && $temp_excl{$group{$t[0]}} && $group{$t[0]} ne "EXCLUDE";
		}
	    }
	}
	$temp_excl{$group{$t[0]}}=$t[2]; # comma separated (per group), to check
	my @excluded;
	if (defined $t[2]){
	    @excluded=split(/\,/,$t[2]);
	}
	foreach my $temp (@excluded){
	    $excl{$group{$t[0]}}{$temp}=1; # real test of the actual excluded tissue
	}
    }
    close GROUPS;
}

open (I, $ARGV[0]);
my $head=<I>;
chomp($head);
my @head=split(/\t/,$head);
my $tis=join("\t",@head[2..$#head]);
my %val; 
for my $i (2..$#head){
    $val{$head[$i]}=$i;
    print "  No group identified for $head[$i]\n" if !$group{$head[$i]} && $groups;
    $group{$head[$i]}=$head[$i] if !$groups;
}

my (%for_matV, %sum, %tot, %max, %min, %OK, %vals);
my (%av, %med, %SD, %min_ev, %max_ev, %TOTAL);
my %tallyTIS_GE;
my %tally;
while (<I>){
    chomp;
    my @t=split(/\t/);
    
    my $g=$t[0];
    my $pre_line=join("\t",@t[0..1]);
    
    for my $i (2..$#t){
	my $tis=$head[$i];
	my $AT;
	if ($group{$tis}){$AT=$group{$tis};} # group
	next if !$AT;
	next if $AT eq "EXCLUDE"; # if the group is EXCLUDE, exclude.
	parse_sample($g,$tis,$AT,@t);
    }
####
    
    ### Parse by tissue groups
    foreach my $tis (sort (keys %{$sum{$g}})){
	if ($tot{$g}{$tis} >= $min_rep){ # i.e. the number of samples for the group is >= min_rep
	    $av{$g}{$tis}=sprintf("%.2f",$sum{$g}{$tis}/$tot{$g}{$tis});
	    $med{$g}{$tis}=sprintf("%.2f",median($g,$tis));
	    $SD{$g}{$tis}=std_dev($g,$tis,$av{$g}{$tis}) if $tot{$g}{$tis}>1;
	    $SD{$g}{$tis}="NA" if $tot{$g}{$tis}==1;
	    $av{$g}{$tis}=$med{$g}{$tis} unless $average; ## default => clean up at some point!!
	    
	    if (defined $min_ev{$g}){
		$min_ev{$g}=$av{$g}{$tis} if $av{$g}{$tis} <= $min_ev{$g};
	    } else {
		$min_ev{$g}=$av{$g}{$tis};
	    }
	    if (defined $max_ev{$g}){
		$max_ev{$g}=$av{$g}{$tis} if $av{$g}{$tis} >= $max_ev{$g};
	    } else {
		$max_ev{$g}=$av{$g}{$tis};
	    }
	    $TOTAL{$g}++;
	}
    }
    
    next if $max_ev{$g} < $min_expr; # filter for minimal expression
    
    foreach my $tis1 (sort (keys %{$av{$g}})){
	my ($real_total, $pos, $neg) = (0,0,0);
	my $sum_others=0;
	my @vals_others=();
	foreach my $tis2 (sort (keys %{$av{$g}})){
	    next if $excl{$tis1}{$tis2}; # 27/03/18
	    
	    if ($tis1 ne $tis2){
		$real_total++; # the total number of tissues matched against tis1
		$sum_others+=$av{$g}{$tis2};
		push(@vals_others,$av{$g}{$tis2});
		if (defined $strict){ # no overlap for any sub-sample
		    $pos++ if (log(($av{$g}{$tis1}+0.01)/($av{$g}{$tis2}+0.01))/log(2) >= log($min_fold)/log(2) 
			       && log(($min{$g}{$tis1}+0.01)/($max{$g}{$tis2}+0.01))/log(2) >= 1
			       && ($min{$g}{$tis1}-$max{$g}{$tis2} >= $min_range));
		    # i.e. log2(2)
		    $neg++ if (log(($av{$g}{$tis1}+0.01)/($av{$g}{$tis2}+0.01))/log(2) <= -1*(log($min_fold)/log(2)) 
			       && log(($max{$g}{$tis1}+0.01)/($min{$g}{$tis2}+0.01))/log(2) <= -1
			       && ($max{$g}{$tis1}-$min{$g}{$tis2} <= -1*$min_range));
		}
		else {
		    $pos++ if (log(($av{$g}{$tis1}+0.01)/($av{$g}{$tis2}+0.01))/log(2) >= log($min_fold)/log(2)) 
			&& ($av{$g}{$tis1}-$av{$g}{$tis2} >= $min_range);
		    $neg++ if (log(($av{$g}{$tis1}+0.01)/($av{$g}{$tis2}+0.01))/log(2) <= -1*(log($min_fold)/log(2))) 
			&& ($av{$g}{$tis1}-$av{$g}{$tis2} <= -1*$min_range);
		}
	    }
	}
	
	$tallyTIS_GE{$tis1}++; # keeps the count of valid tests
	my $dir="";
	$dir="DOWN" if $neg == $real_total; 
	$dir="UP" if $pos == $real_total;
	my $av_others=sprintf("%.2f",$sum_others/$real_total); # the average of the compared tissues
	my $med_others=sprintf("%.2f",median2(@vals_others));
	my $min_others=(sort{$a<=>$b}(@vals_others))[0];
	my $max_others=(sort{$b<=>$a}(@vals_others))[0];
	my $FC_tis1;
	if (defined $average){
	    $FC_tis1=sprintf("%.2f",log(($av{$g}{$tis1}+0.01)/($av_others+0.01))/log(2));
	} else {
	    $FC_tis1=sprintf("%.2f",log(($av{$g}{$tis1}+0.01)/($med_others+0.01))/log(2));
	}
	my $handler="FC_$tis1";
	my $SD_others;
	if ($#vals_others>0){
	    $SD_others=std_dev2(@vals_others) if $#vals_others>0; ### Calculates the SD for "others" (new in V5)
	} else {
	    $SD_others="NA" if $#vals_others==0; # i.e. 1 element
	}
	
	if  (defined $test_tis_OK{$tis1}){
	    if (defined $average){
		print $handler "$pre_line\t$av{$g}{$tis1}\t$av_others\t$FC_tis1\t$SD{$g}{$tis1}\t$SD_others\t$min_others\t$max_others\n";
	    } else {
		print $handler "$pre_line\t$av{$g}{$tis1}\t$med_others\t$FC_tis1\t$SD{$g}{$tis1}\t$SD_others\t$min_others\t$max_others\n";
	    }
	}
	
	if ($dir && abs($FC_tis1)>= log($min_fold_glob)/log(2)){
	    $tally{$tis1}{$dir}++; # keeps the count of tissue-specific
	    if (defined $print){
		if (defined $average){
		    print O "$pre_line\t$tis1\t$dir\t$av{$g}{$tis1}\t$av_others\t$FC_tis1\t$SD{$g}{$tis1}\t$SD_others\t$min_others\t$max_others\n";
		} else {
		    print O "$pre_line\t$tis1\t$dir\t$av{$g}{$tis1}\t$med_others\t$FC_tis1\t$SD{$g}{$tis1}\t$SD_others\t$min_others\t$max_others\n";
		}
	    }
	}
    }
}

print "$v: min_FC >= $min_fold, min_FC_glob >= $min_fold_glob";
print ", min_rep $min_rep";
print ", Min_expr $min_expr";
print ", Min_range $min_range";
print ", Median" if !$average;
print ", Average" if $average;
print "\n";

print "SAMPLE\tTOTAL\tUP\tDOWN\tTot_GE\tPerc_TS-GE\tPerc_UP\tPerc_DOWN\n"; 
my %data_print; my @order;
foreach my $tis (sort (keys %tallyTIS_GE)){ 
    $tally{$tis}{UP}=0 if !$tally{$tis}{UP};
    $tally{$tis}{DOWN}=0 if !$tally{$tis}{DOWN};
    my $t=$tally{$tis}{UP}+$tally{$tis}{DOWN};
    my $fr_GE=sprintf("%.2f",100*$t/$tallyTIS_GE{$tis});
    my $fr_UP=sprintf("%.2f",100*$tally{$tis}{UP}/$tallyTIS_GE{$tis});
    my $fr_DOWN=sprintf("%.2f",100*$tally{$tis}{DOWN}/$tallyTIS_GE{$tis});

    # no longer "-" for downregulated TS
    $data_print{$tis}="$tis\t$t\t$tally{$tis}{UP}\t$tally{$tis}{DOWN}\t$tallyTIS_GE{$tis}\t$fr_GE\t$fr_UP\t-$fr_DOWN\n";
    push(@order,"$fr_GE=$tis");
}
@order = sort {($b=~/(.+?)\=/)[0]<=>($a=~/(.+?)\=/)[0]}(@order);

foreach my $temp (@order){
    my ($tis)=$temp=~/\=(.+)/;
    print $data_print{$tis};
}


###### SUBROUTINES
sub parse_sample {
    my @t2=@_;
    my $g=$t2[0];
    my $tis=$t2[1]; # actual tissuee
    my $AT=$t2[2]; # group
    my @t=@t2[3..$#t2];
    
    # quite convoluted, but maintained
    $sum{$g}{$AT}+=$t[$val{$tis}];
    $tot{$g}{$AT}++;
    
    if (defined $max{$g}{$AT}){
	$max{$g}{$AT}=sprintf("%.2f",$t[$val{$tis}]) if $t[$val{$tis}]>$max{$g}{$AT};
    } else {
	$max{$g}{$AT}=sprintf("%.2f",$t[$val{$tis}]);
    }
    if (defined $min{$g}{$AT}){
	$min{$g}{$AT}=sprintf("%.2f",$t[$val{$tis}]) if $t[$val{$tis}]<$min{$g}{$AT};
    } else {
	$min{$g}{$AT}=sprintf("%.2f",$t[$val{$tis}]);
    }
    $OK{$g}{$AT}.="$tis,";
    $vals{$g}{$AT}.="$t[$val{$tis}],";
}

sub median {
    my @t2=@_;
    my $g=$t2[0];
    my $tis=$t2[1];

    my %V; 
    @{$V{$tis}}=sort{$a<=>$b}(split(/\,/,$vals{$g}{$tis}));
    my $ind="";
    my %med;
    $ind=$#{$V{$tis}}+1; 
    if ($ind%2!=0){
	my $m=int($ind/2);
	$med{$g}{$tis}=$V{$tis}[$m];
    }
    else {
	my $m1=int($ind/2)-1;
	my $m2=int($ind/2);
	$med{$g}{$tis}=($V{$tis}[$m1]+$V{$tis}[$m2])/2;
    }
    return $med{$g}{$tis};
}

sub median2 { # simplified median
    my @vals=@_;
    @vals=sort{$a<=>$b}(@vals);
    my $ind="";
    $ind=$#vals+1;
    my $median2;
    if ($ind%2!=0){
	my $m=int($ind/2);
	$median2=$vals[$m];
    }
    else {
	my $m1=int($ind/2)-1;
	my $m2=int($ind/2);
	$median2=($vals[$m1]+$vals[$m2])/2;
    }
    return $median2;
}

sub std_dev {
    my @t2=@_;
    my $g=$t2[0];
    my $tis=$t2[1];
    my $av=$t2[2];
    if (!$av){
	my $av=sprintf("%.2f",$sum{$g}{$tis}/$tot{$g}{$tis});
    }
    
    my @t3=split(/\,/,$vals{$g}{$tis});
    my $sumCLAS=0;
    my %SD;
    foreach my $n (@t3){
	$sumCLAS+=($n-$av)*($n-$av);
    }
    $SD{$g}{$tis}=sprintf("%.2f",sqrt($sumCLAS/($#t3)));

    return $SD{$g}{$tis};
}

sub std_dev2 {
    my @temp_vals=@_;

    my $sum1=0;
    foreach my $val (@temp_vals){
        $sum1+=$val;
    }
    my $av=$sum1/($#temp_vals+1);
    
    my $sumCLAS=0;
    foreach my $n (@temp_vals){
        $sumCLAS+=($n-$av)*($n-$av);
    }
    my $SD_others=sprintf("%.2f",sqrt($sumCLAS/($#temp_vals)));
    
    return $SD_others;
}
