# $Id$
# FastaReader.pm
# Author Wenchao Lin
# linwenchao@yeah.net
#

package FastaReader;

sub new{
    my $pkg = shift;
    my $file = shift;
    my $headsep = shift;
    my $linesep = shift;

    my $self = {};
    bless $self;

    $self->{headsep} = '>'; 
    $self->{headsep} = $headsep if defined $headsep;
    $self->{linesep} = ''; 
    $self->{linesep} = $linesep if defined $linesep;
    $self->{file} = $file;
    $self->{tell} = CORE::tell($file);

    $self->{buf} = <$file>;
    if (! defined $self->{buf}){
	print STDERR "File appears empty\n";
	return undef;
    }
    if ($self->{buf} !~ /^$self->{headsep}/){
	print STDERR "File doesn't start with a header: $headsep\n";
	return undef;
    }
    chomp $self->{buf};
    return $self;
}

sub getRecord{
    my $self = shift;
    my $head;
    my $data;
    my $file = $self->{file};
    my $tl;

    if (! defined $self->{buf} || $self->{buf} !~ /^$self->{headsep}/){
        return ();
    }

    $head = $self->{buf};
    $head =~ s/^$self->{headsep}//;
    $tl = CORE::tell($file);
    $self->{buf} = <$file>;
    chomp $self->{buf};
    while (defined $self->{buf} && $self->{buf} !~ /^$self->{headsep}/){
        $data .= $self->{buf} . $self->{linesep};
        $tl = CORE::tell($file);
        $self->{buf} = <$file>;
        if (defined $self->{buf}){chomp $self->{buf}};
    }
    $self->{tell} = $tl;
    $self->{seq} = $data;
    $self->{name} = $head;
    return $self;
}


sub seek{
    my $self = shift;
    my $pos = shift;
    my $file = $self->{file};
    my $headsep = $self->{headsep};

    CORE::seek($file, $pos, 0);

    $self->{tell} = CORE::tell($file);
    $self->{buf} = <$file>;
    if (! defined $self->{buf}){
        print STDERR "File appears empty\n";
        return undef;
    }
    if ($self->{buf} !~ /^$self->{headsep}/){
        print STDERR "File doesn't start with a header: $headsep\n";
        return undef;
    }
    chomp $self->{buf};
}


sub tell{
    my $self = shift;    
    return $self->{tell};
}

sub length{
    my $self = shift;
    $self->{len} = length $self->{seq};
    return $self->{len};
}

sub name{
    my $self = shift;
    return $self->{name};
}

sub seq{
    my $self = shift;
    return $self->{seq};
}

1;

=head1 NAME

FastaReader - class for reading a fasta-like formatted file

The assumption is that the file contains a set of multi-line records separated by single-line headers starting with a specific separator (by default '>')

=head1 SYNOPSIS

use FastaReader;

my $file = shift;

open FILE,"<$file";

my $parser = new FastaReader(\*FILE);

while (my $seq = $parser->getRecord()){

    printf("SeqName:%s\n",$seq->name());

    printf("SeqLength:%s\n",$seq->length());

    printf("SeqSequence:%s\n",$seq->seq());

    ...

}

=head1 DESCRIPTION

This module iterates through a fasta-like file retrieving the records in 
(header, data) pairs. When creating a new parser, the user may specify
both the prefix (as a string) that indicates a header (record separator) line
(default is '>'), and a line separator that is used when concatenating the
lines in the input forming the body of each record.  This separator is useful,
for example, when parsing a .qual file where the quality values are separated
by spaces.

=over

=item $parser = new FastaReader($file, $head, $sep) ;

Creates a new parser object reading from file $file and using a specific
header separator ($header - default '>') and a line separator 
($sep - default '')

=item $seq = $parser->getRecord();

Reads a record into an object of obj. If no more records remain returns undef.

=item $name = $seq->name();

Return name of current record

=item $length = $seq->length();

Return length of current record

=item $posn = $seq->tell();

Return offset of current record

=item $parser->seek(posn);

Resets the parser to a specific location (posn) in the file stream.

=back

=head1 AUTHOR

Wenchao Lin E<lt>linwenchao@yeah.netE<gt>

=cut
