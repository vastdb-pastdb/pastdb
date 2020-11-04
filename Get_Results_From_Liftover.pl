#!/usr/bin/perl
# Original name: GetPercCons_AllEV-v2.pl
# Main script used to get the G-conservation by liftOver
# It needs the FILTERED table for each event type
use strict;
use warnings;

my @EVENTS=("EX","INT","ALTA","ALTD");
my @SPECIES=("Aly","Csa","Aal","Bra");

my %cons_ev; # hash of conserved events
my %tot; my %cons; # tallies
my %done_core;

foreach my $test_ev (@EVENTS){
    foreach my $sp (@SPECIES){
	my $in_file="$test_ev-Ath-to-$sp-FILTERED.tab";
	open (LIFT, $in_file) || die "It cannot find liftover output file\n";
	while (<LIFT>){
	    chomp;
	    my @t=split(/\t/);
	    my $ev=$t[0];
	    $cons_ev{$sp}{$ev}=1;
	}
	close LIFT;
	
	open (I, "AllEvents_for_comparison-Ath.txt") || die "Events for comparison";
	<I>;
	while (<I>){
	    chomp;
	    my @t=split(/\t/);
	    my $ev=$t[1];
	    my $reg_type=$t[2];
	    $reg_type="SuperAS" if $t[2] eq "TrulyAS";

	    if ($ev=~/$test_ev/){
		$tot{$test_ev}{$reg_type}{$sp}++;
		$cons{$test_ev}{$reg_type}{$sp}++ if $cons_ev{$sp}{$ev};
		
		if ($reg_type eq "ABI" || $reg_type eq "BIO" || $reg_type eq "TIS"){
		    my $new_type="CORE";
		    unless ($done_core{$sp}{$ev}){
			$tot{$test_ev}{$new_type}{$sp}++;
			$cons{$test_ev}{$new_type}{$sp}++ if $cons_ev{$sp}{$ev};
		    }
		    $done_core{$sp}{$ev}=1;
		}
	    }    
	}
	close I;
    }
}


my @TYPES=("HIGH_PSI","SuperAS","LOW_PSI","CORE","ABI","BIO","TIS");

print "EV_TYPE\tREG_TYPE\tTOTAL_EV";
foreach my $sp (@SPECIES){
    print "\tCONS_$sp\tNO_CONS_$sp\tPercCons_$sp";
}
print "\n";

foreach my $ev_type (@EVENTS){
    foreach my $reg_type (@TYPES){
	print "$ev_type\t$reg_type\t$tot{$ev_type}{$reg_type}{Csa}";
	foreach my $sp (@SPECIES){
	    my $no_cons = $tot{$ev_type}{$reg_type}{$sp} - $cons{$ev_type}{$reg_type}{$sp};
	    my $perc=sprintf("%.2f",100*$cons{$ev_type}{$reg_type}{$sp}/$tot{$ev_type}{$reg_type}{$sp});
	    print "\t$cons{$ev_type}{$reg_type}{$sp}\t$no_cons\t$perc";
	}
	print "\n";
    }
    print "\n";
}


### to get the Fisher CORE vs SuperAS
foreach my $ev_type (@EVENTS){
    foreach my $sp (@SPECIES){
	my $cons_CORE = $cons{$ev_type}{CORE}{$sp};
	my $no_cons_CORE = $tot{$ev_type}{CORE}{$sp} - $cons{$ev_type}{CORE}{$sp};
	my $cons_SAS = $cons{$ev_type}{SuperAS}{$sp};
	my $no_cons_SAS = $tot{$ev_type}{SuperAS}{$sp} - $cons{$ev_type}{SuperAS}{$sp};

	open (TEMP, ">temp.R");
	print TEMP "
data <- matrix(c($cons_CORE, $no_cons_CORE,  $cons_SAS, $no_cons_SAS), nrow = 2, ncol=2, byrow=TRUE)
fisher.test(data,alternative = c(\"two.sided\"))\$p.value\n";
	close TEMP;

	my $p_val=`Rscript ./temp.R`;
	$p_val=~s/\[1\] //g;
	chomp($p_val);
	$p_val=sprintf("%.6f",$p_val);
	
	my $log_ratio = sprintf("%.2f", log(($cons_CORE/$tot{$ev_type}{CORE}{$sp})/($cons_SAS/$tot{$ev_type}{SuperAS}{$sp}))/log(2));
	
	print "$ev_type\t$sp\t$log_ratio\t$p_val\n";
	
	system "rm temp.R";
    }
    print "\n";
}


