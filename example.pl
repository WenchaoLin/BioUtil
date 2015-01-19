#!/usr/bin/perl

use FastaReader;

my $file = shift;

open(IN,$file) || die "Cannot open $ARGV[0]: $!\n";

my $fa = new FastaReader(\*IN);

while (my $seqobj = $fa->getRecord()) {
    printf ("%s\n", $seqobj->seq());
    printf ("%s\n", $seqobj->name());
    printf ("%s\n", $seqobj->length());
}
