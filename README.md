# NAME

FastaReader - class for reading a fasta-like formatted file

The assumption is that the file contains a set of multi-line records separated by single-line headers starting with a specific separator (by default '>')

# SYNOPSIS

use FastaReader;

my $file = shift;

open FILE,"<$file";

my $parser = new FastaReader(\\\*FILE);

while (my $seq = $parser->getRecord()){

    printf("SeqName:%s\n",$seq->name());

    printf("SeqLength:%s\n",$seq->length());

    printf("SeqSequence:%s\n",$seq->seq());

    ...

}

# DESCRIPTION

This module iterates through a fasta-like file retrieving the records in 
(header, data) pairs. When creating a new parser, the user may specify
both the prefix (as a string) that indicates a header (record separator) line
(default is '>'), and a line separator that is used when concatenating the
lines in the input forming the body of each record.  This separator is useful,
for example, when parsing a .qual file where the quality values are separated
by spaces.

- $parser = new FastaReader($file, $head, $sep) ;

    Creates a new parser object reading from file $file and using a specific
    header separator ($header - default '>') and a line separator 
    ($sep - default '')

- $seq = $parser->getRecord();

    Reads a record into an object of obj. If no more records remain returns undef.

- $name = $seq->name();

    Return name of current record

- $length = $seq->length();

    Return length of current record

- $posn = $seq->tell();

    Return offset of current record

- $parser->seek(posn);

    Resets the parser to a specific location (posn) in the file stream.

# AUTHOR

Wenchao Lin <linwenchao@yeah.net>