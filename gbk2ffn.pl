#!/usr/bin/perl
#Convert gbk file into ncbi ffn format 
#genome cds fastas with verbose header containing gi,accession number, genomic start and stop, description
#no copyright Wenchao Lin linwenchao@yeah.net
#Tianjin Biochip Inc.

use warnings;
use diagnostics;
use strict;
use Bio::SeqIO;
use Text::Wrap;

#Fasta file format settings
$Text::Wrap::columns = 80;
$Text::Wrap::separator = "\n";

my $gi_tag = "locus_tag";

if (@ARGV == 0) { die "\tConvert gbk file into ncbi ffn format\n
\tgenome cds fastas with verbose header containing gi,accession number, genomic start and stop, description\n
\tUSAGE: gbk2ffn.pl <gbk> > file.ffn\n
"}

my $seqs = Bio::SeqIO -> new( -file => $ARGV[0], -format => 'genbank');

while (my $seq = $seqs->next_seq) {
    my $accession = $seq -> display_id;
    my $description = $seq->desc;
    for my $feat ($seq -> get_SeqFeatures) {
        if ($feat->primary_tag eq "CDS") {
            my $start = $feat->start;
            my $end = $feat->end;
            #header
            print "\>gi\|";
            if ( $feat->has_tag($gi_tag) ) {
                my @xref = $feat->get_tag_values($gi_tag);
                my $gi = $xref[0];
                $gi =~ s/GI://g;
                print "$gi\|ref\|$accession\|:$start-$end $description";                
            }else{
                print "NA\|ref\|$accession\|:$start-$end $description";
            }
            #sequence:
            print "\n",wrap("","",$feat->seq->seq),"\n";
        }
    }
}
