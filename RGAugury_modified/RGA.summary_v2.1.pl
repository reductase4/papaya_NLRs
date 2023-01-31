#!/usr/bin/perl -w
use strict;
use Scalar::Util qw(looks_like_number);
#count the RGA type and number from RGA.candidates.fas

my %count = ();
my %gene = ();
my $total = 0;

open(IN,$ARGV[0]); # generally the input.RGA.candidates.fasta
while (<IN>) {
    chomp;
    next unless (/^>/);
    my ($gene,$type) = $_ =~ /(.*)\|(.*)/;
    #print "$gene $type\n";
    $count{$type}++;
    $total++;
}
close IN;

my $i = 0;
open(IN,$ARGV[1]); # input.cc_lrr_only.domains.txt
while (<IN>) {
    chomp;
    $i++;
    next if ($i ==1); # ignore the first line which is usually the title of each column.
    my ($id,$len,$lrr,$cc) = split/\t/,$_;
    if ($lrr =~/domain/i and $cc =~/domain/i) {
        my $type = "CL";
        $count{$type}++;
    }
    elsif ($lrr =~/domain/i) {
        my $type = "LRR";
        $count{$type}++;
    }
    elsif ($cc =~/domain/i) {
        my $type = "CC";
        $count{$type}++;
    }
}
close IN;

my @array  = keys %count;

my $printed_count = 0;
#print join("\t","NBS","CNL","TNL","RNL","RN","RX","CN","TN","NL","TX","Others","RLP","RLK","TM-CC\n");
print join("\t","NBS","TX","RX","LRR","CC","CL","CN","TN","RN","NL","TC","RC","TCN","RCN","CNL","TNL","RNL","TCNL","RCNL","RTCNL","Others","RLP","RLK","TM-CC\n");
#foreach my $type ("NBS","CNL","TNL","RNL","RN","RX","CN","TN","NL","TX","OTHER","RLP","RLK","TM-CC") {
foreach my $type ("NBS","TX","RX","LRR","CC","CL","CN","TN","RN","NL","TC","RC","TCN","RCN","CNL","TNL","RNL","TCNL","RCNL","RTCNL","Others","RLP","RLK","TM-CC") {
    if (looks_like_number($count{$type})) {
        print "$count{$type}\t";
        $printed_count += $count{$type};
    }
    else {
        print "0\t";
    }
}
print "\n";
#in case there is other type that is not thoughtfully considered. thus an comparison is set up as ;
my $notication = ($total == $printed_count) ? "\n=>looks good<=\n"  : "\n!!!!!!!!!!!! The total and printed number is unequal !!!!!!!!!!!!\n";
print "$notication\nRGA total = $total;\tPrinted_count = $printed_count\n";
