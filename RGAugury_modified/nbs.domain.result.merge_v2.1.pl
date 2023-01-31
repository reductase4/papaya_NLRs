#!/usr/bin/perl
use strict;
use warnings;
use Getopt::Long;
use Scalar::Util qw(looks_like_number);

GetOptions(my $options = {},
              "-seq=s","-nbs=s","-lrr=s","-tir=s","-rpw8=s","-cc=s","-pfx=s"
);


my $USAGE = <<USAGE;
-------------------------------------------------------------------------------------------------------------
Script: to merge the prediction file

Usage:   perl $0 <options> <files>
-------------------------------------------------------------------------------------------------------------
 coded by Pingchuan Li
Arguments:
        -seq       the fasta seq, to acquire seq length
        -nbs       nbs filename
        -lrr       lrr filename
        -tir       tir filename
        -rpw8      rpw8 filename
        -cc        cc filename  
        -pfx       prefix

enjoy it!
USAGE

my %domain = ();

die $USAGE unless (defined $options->{nbs});

my $fasta   = $options->{seq};
my $nbsfile = $options->{nbs};
my $lrrfile = $options->{lrr};
my $tirfile = $options->{tir};
my $rpw8file = $options->{rpw8};
my $ccfile  = $options->{cc};
my $prefix = ($options->{pfx}) ? $options->{pfx}: "default.";
my $cc_lrr_outfile = $prefix."cc_lrr_only.domains.txt";

domain_sort($nbsfile, "nbs") if (-e $nbsfile);
domain_sort($lrrfile, "lrr") if (-e $lrrfile);
domain_sort($tirfile, "tir") if (-e $tirfile);
domain_sort($rpw8file, "rpw8") if (-e $rpw8file);
domain_sort($ccfile,  "cc" ) if (-e $ccfile );

my %seqlen = ();
open(FASTA,"$fasta");
local $/ = ">";
while (<FASTA>) {
    chomp;
    my ($title,$seq) = split/\n/,$_,2;
    next unless ($title and $seq);
    my ($id) = $title =~ /(\S+)/;
    $seq =~ s/\s+//g;
    my $len = length($seq);
    $seqlen{$id} = $len;
}
close FASTA;

local $/ = "\n";

print "id\tLen\tNBS\tLRR\tTIR\tCC\tRPW8\n";
open(CC_LRR_only, ">$cc_lrr_outfile");
print CC_LRR_only "id\tLen\tLRR\tCC\n";
foreach my $gene (sort {$a cmp $b} keys %domain) {
    my $nbs = ($domain{$gene}->{nbs} and $domain{$gene}->{nbs} =~ /domain/) ? $domain{$gene}->{nbs} : '.' ;
    my $lrr = ($domain{$gene}->{lrr} and $domain{$gene}->{lrr} =~ /domain/) ? $domain{$gene}->{lrr} : '.' ;
    my $tir = ($domain{$gene}->{tir} and $domain{$gene}->{tir} =~ /domain/) ? $domain{$gene}->{tir} : '.' ;
    my $cc  = ($domain{$gene}->{cc}  and $domain{$gene}->{cc}  =~ /domain/) ? $domain{$gene}->{cc}  : '.' ;
    my $rpw8  = ($domain{$gene}->{rpw8}  and $domain{$gene}->{rpw8}  =~ /domain/) ? $domain{$gene}->{rpw8}  : '.' ;

    #make sure all identified genes go through below process as NBS-encoding genes
    #next unless ($nbs =~ /domain/i or $tir =~ /domain/i or $rpw8 =~ /domain/i);
    if ($nbs =~ /domain/i or $tir =~ /domain/i or $rpw8 =~ /domain/i){
        print join("\t", $gene,$seqlen{$gene},$nbs,$lrr,$tir,$cc,$rpw8);
        print "\n";
    }
    elsif ($lrr =~ /domain/i or $cc =~ /domain/i) {
        print CC_LRR_only join("\t",$gene,$seqlen{$gene},$lrr,$cc,"\n");
    }
    
}
close CC_LRR_only;
# -----------subfunctions-----------
sub domain_sort {
    my ($filename,$type) = @_;
    open(IN,$filename);
    while (<IN>) {
        chomp;
        my ($id,@content) = split/\s+/,$_;
        $domain{$id}->{$type} = join(" ",@content);
    }
    close IN;
}
