#!/usr/bin/perl
#extract specific gene sequence from gbk file into ncbi ffn format 
#in-house use for primer design
#no copyright Wenchao Lin linwenchao@yeah.net
##Tianjin Biochip Inc.
#

use warnings;
use diagnostics;
use strict;
use Bio::SeqIO;
use Text::Wrap;

#Fasta file format settings
$Text::Wrap::columns = 80;
$Text::Wrap::separator = "\n";

if (@ARGV == 0) { die "\textract specific gene sequence from gbk file into ncbi ffn format\n\n
\tUSAGE: extract_gene_from_gbk_to_ffn.pl <gbk>  <genename> >gene.ffn\n "}

my $gene_to_extract = $ARGV[1];

print STDERR "reading gbk...\n";
my $seqs = Bio::SeqIO -> new( -file => $ARGV[0], -format => 'genbank');

my $seqid = 0;
while (my $seq = $seqs->next_seq) {
$seqid ++;
print STDERR "processing seq $seqid\r";
    my $accession = $seq -> display_id;
    my $description = $seq->desc;
    for my $feat ($seq -> get_SeqFeatures) {
        if ($feat->primary_tag eq "CDS") {
            my $start = $feat->start;
            my $end = $feat->end;
            if ( $feat->has_tag('gene') ) {
                my @xref = $feat->get_tag_values('gene');
                my $gene = $xref[0];
        if ($gene eq $gene_to_extract){
                print ">${gene}_${accession} $start-$end $description";                
                print "\n",wrap("","",$feat->seq->seq),"\n";
        last;
        }
        }
    }
}
}

