#!/usr/bin/perl

use FastaReader;

my $file = shift;

my $fa = new FastaReader($file);

while (my $seqobj = $fa->next_seq()) {
    printf ("%s\n", $seqobj->seq());
    printf ("%s\n", $seqobj->name);
    printf ("%s\n", $seqobj->length());
}
