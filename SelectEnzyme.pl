#!/usr/bin/perl -w

## Usage:
## Enzyme selection of sequence data
## $0 fasta_file

use strict;
use warnings;

use Getopt::Long;
use List::Util qw(sum);

my $outfile;
my $gapsize;
my $verbose = 1;
my $sizefilter=100;

GetOptions(
  "filter=n" => \$sizefilter,
);


@ARGV == 1 or die "Usage:\n      SelectEnzyme.pl  *.fasta
";
 
my ($fastafile) = @ARGV;

if( defined $gapsize )
{
	if($gapsize !~ m/^\+?\d+\.?\d*$/)
  {
       die "gap size is not a positive size\n";
  }
}

my $make_gaps = 1;
my $begin = 5;
my $end = 20;

my @enzlist = ("NtBspQI","NbBbvCI","NbBsmI","NbBsrDI");
#my @enzlist = ("AflII","BamHI","KpnI","NcoI","NheI","SpeI","BglII","EcoRI","MluI","NdeI","PvuII","XbaI","XhoI");
#my @enzlist = ("AflII","BamHI","KpnI","NcoI","NheI","ApaLI","BglII","EcoRI","HindIII","MluI","NdeI","NotI","PvuII","SpeI","XbaI","XhoI");

my $enzyme_cut_site;
my $cut_site_reg;
my $enz_name ;

my @stats=();
	

print "Enzyme\t5kb-20kb percent(%)\t6kb-15kb percent(%)\t6kb-12kb percent(%)\tAFS(kb)\t#frags>100kb\tMax_frag_size(kb)\n" if($verbose);
print "=========================================================================================================================\n" if ($verbose);

for( my $i=0; $i<= $#enzlist;$i++ )
{
    $enz_name = $enzlist[$i];
    $enzyme_cut_site = _make_enzyme($enz_name);#including cut information
    $cut_site_reg = IUB_to_regexp($enzyme_cut_site); #without cut information
    process_fasta($fastafile,$outfile,$begin,$end,\@stats);
    
#    #print $stats[0]."\n";
   
}
  
  print "=============================Sequence Stats=============================\n";
# print "TotalRealSeqs:". $stats[0]. "\tTotalWithoutNSeqs:".$stats[1];
# printf "\tN Percent:%.2f%%\n",$stats[2];
  #print "TotalSeqs:". $stats[3]. "\tTotalSeqAfterCut:". $stats[4];
  #printf "\tPct:%.2f%%\n", $stats[5];
  #print "Max N gap is ".$stats[6]."\n";
  
  
  fasta_stats($fastafile,$stats[0]/1000000000);


#
#read fasta file
#Streamline processing sequence file which means processing one sequence at one time to remove possible memory problems.
#

sub process_fasta{
	my ($fasta_file,$outfile,$begin,$end,$refstat) = @_;
	
	#my @stats = @{$refstat};
	
	open (FILE,"< $fasta_file") or die("Cannot open file!\n");
	#my %seqs =();


  my $totalseq = 0;
  my $totalrealseq = 0;
  my $totalseqbeforecut = 0;
  my $totalwithoutN = 0;
  
	my $usableseq = 0;
	my $usableseq1 = 0;
	my $usableseq2 = 0;
	my $count100K=0;
	my $count200K=0;
	
	my $header = "";
	my $first = 0;
	my $index = 0;
	#$seqs{$header}="";
	my $sequence;
	my $sequence_withoutN;
	my @afsarray;
	my $max_frag_size=0;
	my $totalsequence = 0;
	my $totalsequence_test = 0;
	my $totalcuts = 0;
	
	my $countN=0;
	my $countn=0;
	
	my @Ngap;
	my $Ngapseq;
	my $maxgap=0;
	my $mingap=0;
	
	$outfile ||= "$fastafile.maps";
  my $fh;
  open $fh,"> $outfile" or die "Can't open '$outfile':$!\n";
  
	my $linecount = 0;
	my $line;
	
  #read fasta sequence file
  
	while( $line = <FILE> ){
		
		 $line =~ s/\r\n$/\n/;
     $line = uc $line;
		chomp($line);
		if ($line =~ /^>/ ){
			$linecount ++;
			if( $linecount > 1)  ##not first line
			{
				 
				 $totalrealseq += length $sequence;
				 
				 @Ngap = $sequence =~ /[Nn]+/g;
				 
				 #$countN += $sequence =~ s/N//g;
				 #$countn += $sequence =~ s/n//g;
				 
	       foreach $Ngapseq (@Ngap)
	       {
	       	  #print $Ngapseq."\n";
				 	  if( length($Ngapseq) > $maxgap )
				    {
				    	 $maxgap = length($Ngapseq);
				    }
				    if( length($Ngapseq) < $mingap and length($Ngapseq) > 1)
				    {
				    	$mingap = length($Ngapseq);
				    }
        }
        
        $sequence_withoutN = $sequence;
				$sequence_withoutN =~ s/[Nn]+//g;
				 
				 
        $totalwithoutN += length $sequence_withoutN;
         my $test1 = length $sequence_withoutN;
         my $test2 = length $sequence;
         
				 
				 $gapsize = 1; # split if N gap is large than 1kb, merge if N gap is less than 1kb
				 my @maps = makeMap($header, $sequence,$enz_name);
				 $totalsequence += (length $sequence_withoutN)/1000;
				 for(my $mapindex =0; $mapindex < @maps; $mapindex++)
				 {
				     if( defined $fh){
				 	      $maps[$mapindex]->write($fh);
				 	      
	  	
	  	   
				 	      my @testfragments = @{$maps[$mapindex]->fragments};
				 	      
				 	      for(my $k=0; $k <= $#testfragments;$k++ )
	  	          {
	  		          if( $testfragments[$k] > 100 )
				          {
				 	      	   	$count100K ++;
				          }
				 	      	
				          if( $testfragments[$k] > 200 )
				          {
				 	      	  	$count200K ++;
				          }
				 	      	
				          if( $testfragments[$k] > $max_frag_size )
				          {
				 	      	  	$max_frag_size = $testfragments[$k];
				          }
	
	  	          }
	
				 	      
				 	      $totalseqbeforecut += sum(@testfragments);
				 	     
				 	     # $totalsequence += sum(@testfragments);
				        $totalcuts += $#testfragments;
				 	      
				 	      if( sum(@testfragments) > $sizefilter and $#testfragments > 0  )
				 	      {
	                			 	        
				 	        $totalseq += sum(@testfragments);
				 	        my $afs;
				 	      
				 	        $afs = sum(@testfragments)/($#testfragments);
				 	       #$totalsequence += (length $sequence)/1000;
				 	       #$totalcuts += @testfragments;
				 	       
				 	       
				 	       
				 	       #$afs = (length $sequence)/(1000*(@testfragments));
				 	       
				 	        if( $afs >= $begin and $afs <= $end )
				 	        {
				 	      	  $usableseq += sum(@testfragments);
				 	        }
				 	        if( $afs >= 6 and $afs <= 12 )
				 	        {
				 	      	  $usableseq1 += sum(@testfragments);
				 	        }
				 	        if( $afs >=6 and $afs <= 15 )
				 	        {
				 	      	  $usableseq2 += sum(@testfragments);
				 	        }
				 	      
				 	        push @afsarray, $afs;
				 	      }
				     }
				 }
				 @maps = ();
				 
				 #print $totalwithoutN."\n";
				 
		     $index++;
			}
			
			#decriptive header in fasta file
		  
			$header = $line;
			$header =~ s/^>//;
      $header =~ s/\s.*//g;
      $header =~ s/(:|\&|\(|\)|,|\\|\/|\#|\!)/\_/g;  #processing some strange characters :&{},\/#
			
			#print $header."\n";
			
			$sequence = "";#important to remove uninitialized warning.
			if ($first == 0){
				$first = 1;
			}
			next;
		}
		if($first == 0){
			die ("NOT a standard fasta file.\n");
		}
		else
		{ 
	        $sequence=$sequence.$line; 
	  }
	}
	
	$totalrealseq += length $sequence;
	
	 $countN += $sequence =~ s/N//g;
	 $countn += $sequence =~ s/n//g;
	
	  foreach $Ngapseq (@Ngap)
	   {
				 	  if( length($Ngapseq) > $maxgap )
				    {
				    	 $maxgap = length($Ngapseq);
				    }
				    if( length($Ngapseq) < $mingap and length($Ngapseq) > 1)
				    {
				    	$mingap = length($Ngapseq);
				    }
     }
	
	  $sequence_withoutN = $sequence;
	  $sequence_withoutN =~ s/[Nn]+//g;
	  
	  $totalwithoutN += length $sequence_withoutN;
	  
 
	$gapsize = 1;
	my @maps = makeMap($header, $sequence,$enz_name);
	$totalsequence += (length $sequence_withoutN)/1000;
	
	for(my $mapindex =0; $mapindex < @maps; $mapindex++)
	{
				 $maps[$mapindex]->write($fh);
				 my @testfragments = @{$maps[$mapindex]->fragments};
				 
				 for(my $k=0; $k <= $#testfragments;$k++ )
	  	   {
	  		          if( $testfragments[$k] > 100 )
				          {
				 	      	   	$count100K ++;
				          }
				 	      	
				          if( $testfragments[$k] > 200 )
				          {
				 	      	  	$count200K ++;
				          }
				 	      	
				          if( $testfragments[$k] > $max_frag_size )
				          {
				 	      	  	$max_frag_size = $testfragments[$k];
				          }
	
	  	  }

				 
				 #$totalseqbeforecut += sum(@testfragments);
				 
				 #$totalsequence += sum(@testfragments);
				 $totalcuts += $#testfragments;
				 	  
				 
				 if( sum(@testfragments) > $sizefilter and $#testfragments > 0 )
				 {
				 	  $totalseq += sum(@testfragments);
				  
				    my $afs;
				   $afs = sum(@testfragments)/$#testfragments;
				   
				 
				    if( $afs >= $begin and $afs <= $end )
				    {
				 	      	$usableseq += sum(@testfragments);
				    }
				    if( $afs >= 6 and $afs <= 12 )
				    {
				 	      	$usableseq1 += sum(@testfragments);
				    }
				    if( $afs >=6 and $afs <= 15 )
				    {
				 	      	$usableseq2 += sum(@testfragments);
				    }
				 
				    push @afsarray, $afs;
				  }
				 
	}
	
	#print "Total sequence is:".$totalsequence." ".$totalcuts."\n";		 			 	 
  my $afsnew;
  if( $totalcuts > 0 )
  {
  	$afsnew = $totalsequence/$totalcuts;
  }
  
#print $countN." ".$countn."\n";
if( $totalseq > 0 )  
{
  	printf "%s\t%13.2f\t%17.2f\t%21.2f\t%13.2f\t%7d\t%15.2f\n",$enz_name,$usableseq/$totalseq*100,  $usableseq2/$totalseq*100,$usableseq1/$totalseq*100,
  $afsnew,$count100K, $max_frag_size;
}
else
{
	  printf "%s\tMax contig size without N is less than 100kb\n",$enz_name;
}
	
  $refstat->[0] = $totalrealseq;
  $refstat->[1] = $totalwithoutN;
  $refstat->[2] = ($totalrealseq-$totalwithoutN)/$totalrealseq * 100;
  $refstat->[3] = $totalseqbeforecut;
  $refstat->[4] = $totalseq;
  if( $totalseqbeforecut > 0 )
  {
  	$refstat->[5] = $totalseq/$totalseqbeforecut*100;
  }
  $refstat->[6] = $maxgap;
  
  #@stats = ($totalrealseq,$totalwithoutN,($totalrealseq-$totalwithoutN)/$totalrealseq * 100, $totalseqbeforecut, $totalseq, $totalseq/$totalseqbeforecut*100);
  
  #print "TotalRealSeqs:". $totalrealseq. " TotalWithoutNSeqs:".$totalwithoutN;
  #printf " N Percent:%.2f%%\n",($totalrealseq-$totalwithoutN)/$totalrealseq * 100;
  #print "TotalSeqs:". $totalseqbeforecut. " TotalSeqAfterCut:". $totalseq;
  #printf " Pct:%.2f%%\n", $totalseq/$totalseqbeforecut*100;
 
  
	$index++;
	
	$outfile =~ s/ /\\ /g;
	system("rm $outfile");
	
	#print "Total number of maps is " . $index. "\n";
	close(FILE);  
}

#
#Read fasta file
#Return sequence reference of hash variable

sub read_fasta{
	open (FILE,"<@_") or die("Cannot open @_:$!\n");
	my %seqs =();
	my $header = "";
	my $first = 0;
	my $Ncount=0;
	my $Acount=0;
	my $Ccount=0;
	my $Tcount=0;
	my $Gcount=0;
	my $Unknowncount=0;
	my @NAcount;
	my $TotalNA=0;
	my $CharCount;
	my $index = 0;
	my $sequence;
	my $totalcontigCount=0;
	my $maxcontigLength=0;
	my $totalcontigsize =0;
	my $numberofsequences =0;
	
	#$seqs{$header}=();
	
	my $fh;
	if( defined $outfile)
	{
		open $fh,"> $outfile" or die "Can't open '$outfile':$!\n";
	}
	
	my $linecount = 0;
	
	while( my $line = <FILE>){
		
		$line = uc $line;
		chomp($line); 
		
		if ($line =~ /^>/ ){
			$linecount ++;
			$numberofsequences++;
			if( $linecount > 1)  ##not first line
			{
			
			@NAcount = countNA($sequence);
			$Ncount += $NAcount[0];
			#$Acount += $NAcount[1];
			#$Ccount += $NAcount[2];
			#$Tcount += $NAcount[3];
			#$Gcount += $NAcount[4];
			#$Unknowncount += $NAcount[5];
			$TotalNA += length $sequence ;
			
			if(defined $outfile)
			{
				my $maskoutput = generate_mask($sequence);
				print $fh $header." ".$maskoutput."\n";
			}	
			$CharCount= [$Ncount, $Acount,$Ccount, $Tcount,$Gcount,$Unknowncount,$TotalNA];
			$seqs{$header}= $CharCount;
				
			$index ++;
			}
			
			
			#decriptive header in fasta file
			$header = $line;
			$header =~ s/^>//;
			$header =~ s/\s.*//;
			
				
			$Ncount = $Acount = $Ccount = $Tcount= $Gcount = $Unknowncount= 0;
			$TotalNA = 0;
			$sequence = "";
			
      if ($first == 0){
				$first = 1;
			}
			next;
		}
		if($first == 0){
			die ("Not a standard Fasta file.\n");
		}
		else{	
      $sequence=$sequence.$line; 
    }
	}

  #last sequence
        
  @NAcount = countNA($sequence);
	$Ncount += $NAcount[0];
	$TotalNA += length $sequence ;
			
  $CharCount= [$Ncount, $Acount, $Ccount, $Tcount, $Gcount, $Unknowncount, $TotalNA];
	$seqs{$header}=$CharCount; # last line
        
  # print $header." ".$TotalNA."\n";	
  printf "%23s\t%7d\n","Number of contigs:",$numberofsequences;
  
	close(FILE);
	return \%seqs;
}


#
#Generate fasta stats
#

sub fasta_stats{
	my ($fastafile,$genomesize) = @_;
	
	my $seqs = &read_fasta($fastafile);
	

  my @array =  keys %$seqs;
 # print "Number of scaffold sequences:".@array."\n";

  my $TotalNA = 0;
  my $TotalN = 0;

  my $PrintN50 = 0;
  my $PrintN80 = 0;
  my $PrintN90 = 0;
  my $count = 0;

  my $size200K= 0;
  my $size100K = 0;
  my $size50K = 0;
  my $size20K = 0;
  
  my $total200K;
  my $count200K;

  my @sortarray = sort {$$seqs{$b}->[6] <=> $$seqs{$a}->[6]} keys %$seqs;

  foreach(sort {$$seqs{$b}->[6] <=> $$seqs{$a}->[6]} keys %$seqs){       #descend sorting based on total NA
    $TotalN += $$seqs{$_}->[0];
    $TotalNA += $$seqs{$_}->[6];
    $count++;
    
	if( $TotalNA >= $genomesize *1000000000 * 0.5 && $PrintN50 == 0 ){
	   printf "%23s\t%7.2fkb\n","N50:",$$seqs{$_}->[6]/1000;
	   $PrintN50 = 1;
	}
	if( $TotalNA >= $genomesize *1000000000 * 0.9 && $PrintN80 == 0){
	   printf "%23s\t%7.2fkb\n","N90:",$$seqs{$_}->[6]/1000;
	   $PrintN80 = 1;
	}
	if( $TotalNA >= $genomesize *1000000000 * 0.95 && $PrintN90 == 0 ){
	   printf "%23s\t%7.2fkb\n","N95:",$$seqs{$_}->[6]/1000;
	   $PrintN90 = 1;
	}

        my $printcount = $count -1;
        
        if( $$seqs{$_}->[6] < 200*1000 and $size200K == 0)
        {
           $total200K = $TotalNA-$$seqs{$_}->[6];
           $count200K = $printcount;
           #print "number of scaffold larger than 200K is ".$printcount." ".$TotalNA."\n";
           $size200K = 1;
        }
        if( $$seqs{$_}->[6] < 100* 1000 and $size100K == 0)
        {
          # print "number of scaffold larger than 100K is ".$printcount." ".$TotalNA."\n";
           $size100K = 1;
        }
        if( $$seqs{$_}->[6] < 50*1000 and $size50K == 0 )
        {
          #  print "number of scaffold larger than 50K is ".$printcount." ".$TotalNA."\n";
            $size50K = 1;
        }
        if( $$seqs{$_}->[6] < 20*1000 and $size20K == 0 )
        {
          #  print "number of scaffold larger than 20K is ".$printcount." ".$TotalNA."\n";
            $size20K = 1;
        }
        

 }

  if( $size200K == 0 )
  {
        	 $total200K=$TotalNA;
        	 $count200K=$count;
        	 $size200K=1;
  }

printf "%23s\t%7.2fkb\n","Min contig size:",$$seqs{$sortarray[$#sortarray]}->[6]/1000;
printf "%23s\t%7.2fkb\n","Max contig size:",$$seqs{$sortarray[0]}->[6]/1000;
printf "%23s\t%7.2fkb\n","Ave contig size:",($TotalNA/$count)/1000;
printf "%23s\t%7.2f%%\n","Contig>200kb percent:",$total200K/$TotalNA*100;
printf "%23s\t%7d\n","Contig>200kb count:",$count200K;
printf "%23s\t%7.2f%%\n","N percent:",$TotalN/$TotalNA*100; 
printf "%23s\t%7.2fMb\n","Total size:",$TotalNA/1000000;

}

#
#Count number of NAs in each sequence.
#

sub countNA{
	my ($line) = @_;
	my @count = ();
	my $Acount=0;
	my $Ccount=0;
	my $Tcount=0;
	my $Gcount=0;
	my $Ncount=0;
	my $Unknowncount=0;
	my $Allcount = 0;
	
	chomp($line);
	
	$Ncount = ($line =~ s/N//g);
	#$Acount = ($line =~ s/A//g);
	#$Ccount = ($line =~ s/C//g);
	#$Tcount = ($line =~ s/T//g);
	#$Gcount = ($line =~ s/G//g);
	#$Unknowncount = ($line =~ s/[^ACTGN]//g);
	
=begin
	if($Ncount eq ""){
	 	$Ncount = 0;
	 }
	 
	if($Acount eq ""){
	 	$Acount = 0;
	 }   
	
	if($Ccount eq ""){
	 	$Ccount = 0;
	 }
	 
	 if($Tcount eq ""){
	 	$Tcount = 0;
	 } 
	   
	 if($Gcount eq ""){
	 	$Gcount = 0;
	 }
	 
	 if($Unknowncount eq ""){
	 	$Unknowncount = 0;
	 }
=cut 
	 @count = ($Ncount,$Acount,$Ccount,$Tcount,$Gcount,$Unknowncount);
	   
	return @count;
}

#
#Generate optical maps
#

sub makeMap
{
  my ($header, $seq, $enz,$enzyme_cut_site) = @_;
  
  my $fragmentarray = _seq2fragments($seq,$enz,$enzyme_cut_site);
 
  my @maps;
  my $chunk = 0;
  
  for(my $i= 0; $i < @$fragmentarray; $i++ )
  {
    if( defined $gapsize )
  	{
  		push @maps, Maps->new(
      id        => $header."_".$chunk,
      enzyme    => $enz,
      fragments => $fragmentarray->[$i]
      );
    $chunk++;
    }
    else
    {
    	push @maps, Maps->new(
      id        => $header,
      enzyme    => $enz,
      fragments => $fragmentarray->[$i]
      );
    }
    
  
  }
  return @maps;
}

#
#Enzyme generation without BioPerl
#

sub _make_enzyme
{
  #enzyme database, borrowed from Deacon's enzyme stats
  my %enzymes = (
  AflII   => { cut_site => 'C^TTAAG', },
  BamHI   => { cut_site => 'G^GATCC', },
  KpnI    => { cut_site => 'GGTAC^C', },
  NcoI    => { cut_site => 'C^CATGG', },
  NheI    => { cut_site => 'G^CTAGC', },
  ApaLI   => { cut_site => 'G^TGCAC', },
  BglII   => { cut_site => 'A^GATCT', },
  EcoRI   => { cut_site => 'G^AATTC', },
  HindIII => { cut_site => 'A^AGCTT', },
  MluI    => { cut_site => 'A^CGCGT', },
  NdeI    => { cut_site => 'CA^TATG', },
  NotI    => { cut_site => 'GC^GGCCGC', },
  PvuII   => { cut_site => 'CAG^CTG', },
  SpeI    => { cut_site => 'A^CTAGT', },
  XbaI    => { cut_site => 'T^CTAGA', },
  XhoI    => { cut_site => 'C^TCGAG', },
  AatII   => { cut_site => 'GACGT^C', },
  Acc65I  => { cut_site => 'G^GTACC', },
  AgeI    => { cut_site => 'A^CCGGT', },
  ApaI    => { cut_site => 'GGGCC^C', },
  AseI    => { cut_site => 'AT^TAAT', },
  AvrII   => { cut_site => 'C^CTAGG', },
  BmtI    => { cut_site => 'GCTAG^C', },
  BseYI   => { cut_site => 'C^CCAGC', },
  BspEI   => { cut_site => 'T^CCGGA', },
  BspHI   => { cut_site => 'T^CATGA', },
  BsrGI   => { cut_site => 'T^GTACA', },
  BssSI   => { cut_site => 'C^ACGAG', },
  Bsu36I  => { cut_site => 'CC^TNAGG', },
  EagI    => { cut_site => 'C^GGCCG', },
  KasI    => { cut_site => 'G^GCGCC', },
  MfeI    => { cut_site => 'C^AATTG', },
  NsiI    => { cut_site => 'ATGCA^T', },
  PaeR7I  => { cut_site => 'C^TCGAG', },
  PciI    => { cut_site => 'A^CATGT', },
  PspOMI  => { cut_site => 'G^GGCCC', },
  PstI    => { cut_site => 'CTGCA^G', },
  SacI    => { cut_site => 'GAGCT^C', },
  SalI    => { cut_site => 'G^TCGAC', },
  SphI    => { cut_site => 'GCATG^C', },
  TspMI   => { cut_site => 'C^CCGGG', },
  XmaI    => { cut_site => 'C^CCGGG', },
  NtBspQI => { cut_site => 'GCTCTTCN^N'},
  NbBbvCI => { cut_site => 'CCTCA^GC'},
  NbBsmI  => { cut_site => 'GAATG^CN'},
  NbBsrDI => { cut_site => 'GCAATG^NN'},
);
  
  my $enz = shift;
  my $cutsite;
 
  for my $item (keys(%enzymes) )
  {
    if ($item eq $enz ){
        $cutsite = $enzymes{$item}->{cut_site};
        return $cutsite;
    }
  }
  die "Unknown enzyme: $enz\n";
 }

#
#generate sequence fragments
#

sub _seq2fragments
{
  my ($seq, $enz,$enzyme_cute_site) = @_;
  my $cut_after = index($enzyme_cut_site,'^');
  
  #including cut location in reporting fragment size     
  my ($site_3prime_seq,$site_5prime_seq);   
  if( $cut_after == 0 )
  {
  	$site_3prime_seq ="";
  	$site_5prime_seq =$seq;	
  }
  elsif($cut_after == length $enzyme_cut_site){
  	$site_3prime_seq = $seq;
  	$site_5prime_seq = "";
  }
  else{
  	$site_3prime_seq = substr($enzyme_cut_site,0,$cut_after);
  	$site_5prime_seq = substr($enzyme_cut_site,$cut_after+1,length($enzyme_cut_site)-$cut_after);
  }

  my @seqcuts;
  my @array; 
  if( $make_gaps )
  {  
    my @chunks = split( /(N+)/i, $seq);
    foreach( @chunks )
    {
      if( $_ =~ /N/i ) # don't bother cutting if N
      {
        push @seqcuts, $_;
      }
      elsif ( $_ =~ /./ ) # make sure not null
      {
        @array = split(/$cut_site_reg/,$_, -1);
        
       for( my $i=0; $i < scalar @array; $i++ )
       {
          	$i < $#array and $array[$i]=$array[$i].$site_3prime_seq;
  	        $i > 0 and $array[$i]=$site_5prime_seq.$array[$i];
        }
       push @seqcuts, @array;
       }
    }
  }
  else
  {
  
    @array = split(/$cut_site_reg/,$_, -1);
        
       for( my $i=0; $i < scalar @array; $i++ )
       {
          	$i < $#array and $array[$i]=$array[$i].$site_3prime_seq;
  	        $i > 0 and $array[$i]=$site_5prime_seq.$array[$i];
        }
     push @seqcuts, @array;
  }
  

  my $index;
  my @fragmentarray; #{\@fragments1,\@fragments2,\@fragments3,...}
  my @fragments;     #{fragment1,fragment2,fragment3,...}
  my $li;
  
  if(defined $gapsize) #split seqcuts based on defined value
  {
  	
  	@fragmentarray = _split_filter(\@seqcuts,$gapsize); #this fragmentarray contains sequence information
  	@fragmentarray = _merge_filter(\@fragmentarray,$gapsize); #this fragmentarray contains fragment size information
  	
  }
  else
  {

#  	for($index = 0; $index <= $#seqcuts; $index ++ )
#  	{ 
#  		$li=length($seqcuts[$index])/1000;
#  		$li .= "G" if $seqcuts[$index] =~ /^N+$/i and $make_gaps;
#  		push @fragments,$li;
#  		
#  	}
#  	push @fragmentarray,\@fragments;

   #print $#seqcuts."\n";
   push @fragmentarray,\@seqcuts;
   @fragmentarray = _merge_filter(\@fragmentarray,1000000);   #always merge
   
  }
  
  return \@fragmentarray;
}

#
#standard IUB code from Decon's code
#

sub IUB_to_regexp {
  my ($iub) = @_;
  my $regular_expression = '';
  my %iub2character_class = (
        A => 'A',
        C => 'C',
        G => 'G',
        T => 'T',
        R => '[GA]',
        Y => '[CT]',
        M => '[AC]',
        K => '[GT]',
        S => '[GC]',
        W => '[AT]',
        B => '[CGT]',
        D => '[AGT]',
        H => '[ACT]',
        V => '[ACG]',
        N => '[ACGT]',
    );
  my %reviub2character_class = (
        A => 'T',
        C => 'G',
        G => 'C',
        T => 'A',
        R => '[CT]',
        Y => '[GA]',
        M => '[TG]',
        K => '[CA]',
        S => '[GC]',
        W => '[AT]',
        B => '[CGA]',
        D => '[ACT]',
        H => '[AGT]',
        V => '[TCG]',
        N => '[ACGT]',
    );
  $iub =~ s/\^//g;
  my $re1 = '';
  for (my $i=0; $i<length($iub); ++$i) {
    $re1 .= $iub2character_class{substr($iub, $i, 1)};
  } 
  my $re2 = '';
  for (my $i=length($iub)-1; $i>=0; --$i){
    $re2 .= $reviub2character_class{substr($iub, $i, 1)};
  }
  if ($re1 ne $re2){
    $regular_expression = $re1 . '|' . $re2;
  } else {
    $regular_expression = $re1;
  }
  return $regular_expression;
}


#
#merge fragment based on minimum size and return fragment size array
#

sub _merge_filter
{
  my ($seqcuts, $min) = @_;

  #print $min."\n";
  my $fragments = [@$seqcuts];

  my @fragmentarray;
  
  my $i;
  my $j;
  
  #print @$fragments."\n";
  
  for($i=0; $i< @$fragments; $i++ )
  {
   my $fragmentpiece = $fragments->[$i];
   my $fragmentpiecelength;
   
   for($j=0; $j < @$fragmentpiece; $j++)
   {
   	#    print "debug pos 5:".$j." ".$fragmentpiece->[$j]."\n";
  	    if( $fragmentpiece->[$j] =~ /N/i )
        {
        	 #print "debug pos 6:".$fragmentpiece->[$j]."\n";
             my $size = length($fragmentpiece->[$j])/1000;
             
             next if( $size > $min);  
             splice @$fragmentpiece, $j, 1; #remove the gap
  
             my $index = $j-1;
             my $next;
             my $lengthnext;

            if( $index < 0 )
             { $index = 0;  } 
            else
             { $next = splice @$fragmentpiece, $j, 1;  }

             $next = "" unless defined $next;
             $lengthnext = length($next)/1000;
              
             $fragmentpiecelength->[$index] += $size + $lengthnext;
             $j--; # decrement so that the for loop worksvi
          }
          else
          {
          	my $size = length($fragmentpiece->[$j])/1000;
          	
          	$fragmentpiecelength->[$j] = $size;
          }
     } # end if gap
     push @fragmentarray,$fragmentpiecelength;
  } # end for all fragments
  
  return @fragmentarray;
}

#
#split fragments based on gap size.
#output is fragmentarray with reference to each fragment (sequence information, not size informaiton)
#

sub _split_filter
{
  my ($seqcuts, $min) = @_;

  my $fragments = [ @$seqcuts ];

  my @fragmentarray;
  
  for (my $i=0; $i<= $#$fragments; $i++)
  {
  	
    if( $fragments->[$i] =~ /N/ )
    {
      
      my $size = length($fragments->[$i])/1000;

      next if( $size < $min); 

      splice @$fragments, $i, 1; # remove the gap
      next if $i<0; 

      # take the rest of the data
      my @pre = splice @$fragments, 0, $i;
      $i=-1;
 
      next unless scalar( @pre ) > 0;
      
      my @prelength;
      
      for(my $j=0; $j <= $#pre; $j++)
      {
      	$prelength[$j] = length($pre[$j])/1000;
      }
      
      push @fragmentarray,\@pre;

    } # end if gap
  } # end for all fragments
  
  my @lastlength;
  my @lastpre=@$fragments;
  
  for(my $j=0; $j <= $#$fragments;$j++)
  {
  	$lastlength[$j] = length($lastpre[$j])/1000;
  }

  if(@lastlength)   #last fragment is N fragment 
  {
  	 push @fragmentarray,\@lastpre;
  }
 
  #print "fragment array count is ".@fragmentarray."\n";

  return @fragmentarray;
}

#==================================
# Packages
#==================================

package Maps;

use strict;
use warnings;

=head2 new

Create a new OpGen::MapDB::MapData::Map object.  The constructor
takes a possibly empty hash reference of attributes.

  my $map = OpGen::MapDB::MapData::Map->new(
    id        => 'mapdb:1',
    enzyme    => 'XhoI',
    fragments => [ 1.1, 2.2, 3.3 ]
  );

The 'id' attribute is the map identifier used in the mapset.  The
'enzyme' is the name of the restriction enzyme.  The 'fragments'
attribute is a list of fragments represented in Kbp.

=cut

sub new {
  my $class = shift;
  my $self  = {
    id        => undef,
    enzyme    => undef,
    fragments => undef,
    is_circular => 0,
    @_
  };
  return bless($self, $class);
}


#
# private function to allow easy setting/retrieval of parameters
#
sub _parameter
{
  my ($self, $token, $val) = @_;
  $self->{$token} = $val if defined $val;
  $self->{$token};
}

=head2 id

Set or return the map identifier attribute:

  $map->id('mapdb:1');
  my $id = $map->id();

=cut

sub id {
  my $self = shift;
  $self->_parameter("id", @_);  
}

=head2 type

Set or return the type attribute:

  $map->type('consensus');
  my $type= $map->type();

=cut

sub type {
  shift->_parameter("type", @_);  
}

=head2 enzyme

Set or return the enzyme attribute:

  $map->enzyme('BahHI');
  my $enzyme = $map->enzyme();

=cut

sub enzyme {
  my $self = shift;
  $self->_parameter("enzyme", @_);  
}

=head2 fragments

Set or return the fragments attribute (a list reference):

  $map->fragments([ 1.1, 2.2, 3.3 ]);
  my $fragments = $map->fragments();

=cut

sub fragments {
  my $self = shift;
  $self->_parameter("fragments", @_);  
}

=head2 is_circular

Get/set if this map is circular

=cut
sub is_circular
{
  shift->_parameter('is_circular', @_)
}

=head2 number_fragments

Return the number of fragments, which includes the number of gaps, since
gaps are regarded as special kinds of fragments.

  my $number = $map->number_fragments();

=cut

sub number_fragments {
  my $self = shift;

  my $rvalue = (scalar @{ $self->fragments() } );
  return $rvalue;
}

=head2 number_gaps

Return the number of gaps

  my $number = $map->number_gaps();

=cut

sub number_gaps{
  my $self = shift;

  my $g=0;
  foreach ( @{ $self->fragments() } )
  {
    $g++ if /G/;
  }
  return $g;
}
=head2 size

Return the size of the map, including any gaps

  my $size = $map->size();

=cut

sub size {
  my $self = shift;
  my $s=0;
  for ( @{ $self->fragments() } )
  {
    (my $v = $_ ) =~ s/G//;
    $s += $v;
  }
  return $s;
}

=head2 write

Write map data to the specified file handle.  See text for the form of
the output:

  $map->write($fh);

=cut

sub write {
  my $self = shift;
  my $fh   = shift;

  print $fh $self->text;
}

=head2 text

Return the map data to a variable in the same form as if it were written
to a file with write.  The form of the output is:

<id>
<tab> <enzyme> <tab> <first letter of enzyme> <tab> <tab-separated fragments>
<blank line>

=cut
sub text {
  my $self = shift;

  my $text = $self->id() . "\n";
  $text .= "\t" . $self->enzyme() . "\t" . substr($self->enzyme(), 0, 1);

  foreach my $f (@{$self->fragments()}) {
    if( $f =~ /G/ )
    {
      $f =~ s/G//;
      $text .= sprintf("\t%.03f", $f) . "G";
    }
    else
    {
      $text .= sprintf("\t%.03f", $f);
    }
  }

  $text .= "\n\n";
  return $text;
}

=head2 read

Class method read() reads data from a file handle and return a map object.
If successful, this method returns a map object; otherwise, this method
returns undef.

=cut

sub read {
  my $class = shift;
  my $fh    = shift;

  while (<$fh>) {
    chomp;
    last if !/^(#.*|\s*)$/;
  }

  my $id   = $_;
  my $data = <$fh>;

  return undef if (!defined($id) || !defined($data));

  my @list = split(/\s+/, $data);
  shift @list;
  my $enzyme = shift @list;
  shift @list;
  my @fragments = ( );
  foreach (@list) {
    my ($kb) = /^([\d\.]+).*/;
    $kb .= 'G' if /G/;
    push(@fragments, $kb );
  }

  return undef if (!defined($enzyme) || (scalar(@fragments) == 0));

  return Maps->new(
    id        => $id,
    enzyme    => $enzyme,
    fragments => \@fragments
  );
}

=head2 parameter

Get or set any arbitrary parameter.  This routine allows applications to store
additional meta-data in the map object.

  $map->parameter('note', 'Created by me');
  my $note = $map->parameter('note');

=cut
sub parameter
{
  _parameter(@_);
}

=head2 scale_fragments

Scale all the fragments by the given factor.  Useful for converting from
kilobases to bases and back.

=cut
sub scale_fragments
{
  my ($self, $scale) = @_;

  my $frags = $self->fragments;

  for my $i ( 0 .. $#$frags)
  {
    if ($frags->[$i] =~ /G/)
    {
      (my $a = $frags->[$i] ) =~ s/G//;
      $frags->[$i] = $a * $scale;
      $frags->[$i] .= "G";
    }
    else
    {
      $frags->[$i] *= $scale; 
    }
  }
}

1;

package MapSets;

use strict;
use warnings;

=head2 new

Create a new OpGen::MapDB::MapData::MapSet object.  The constructor
takes no arguments.

  my $mapset = OpGen::MapDB::MapData::MapSet->new();

=cut

sub new {
  my $class = shift;
  my $self  = { maps => {} };
  return bless($self, $class);
}

=head2 add

Add a new map (or a list of maps) to the mapset.  Note that each map
must have a unique identifier.  If two maps share the same identifier,
then only one of the two maps will be stored in the database.

  $mapset->add($map);
  $mapset->add($map1, $map2);
  $mapset->add($mapset2->maps());
  $mapset->add(@maplist);

=cut

sub add {
  my $self = shift;
  my @maps = @_;

  foreach (@maps) {
    $self->{'maps'}{$_->id()} = $_;
  }
}

=head2 remove

Remove the maps (specified by map identifier) from the mapset.
Map identifiers for those maps not in the mapset are ignored.

  $mapset->remove($map->id());
  $mapset->remove($mapset2->ids());

=cut

sub remove {
  my $self   = shift;
  my @remove = @_;

  foreach (@remove) {
    delete $self->{'maps'}{$_} if defined($self->{'maps'}{$_});
  }
}

=head2 map

Return the map associated with the specified map identifier.
If the specified map does not exist in the mapset, then undef
is returned.

  my $map  = $mapset->map($id);

=cut

sub map {
  my $self = shift;
  my $id   = shift;

  return $self->{'maps'}{$id};
}

=head2 maps

Return all of the maps in this mapset in a list of Map objects.  You can
also specify a list of ids, in which case it will only return the matching
ids.

  my @maps = $mapset->maps();       # get all maps
  my @maps = $mapset->maps(@ids);   # get subset

=cut

sub maps {
  my $self = shift;
  my @ids = @_;

  if( @ids ) {
    my @rvalue = ();
    foreach (@ids) {
      push @rvalue, $self->map($_);
    }
    return @rvalue;
  }
  else {
    return values(%{$self->{'maps'}});
  }
}

=head2 ids

Return all of the map ids associated with the maps in this mapset.

  my @ids  = $mapset->ids();

=cut

sub ids {
  my $self = shift;

  return keys(%{$self->{'maps'}});
}

=head2 number_maps

Return the number of maps in the mapset.

=cut
sub number_maps {
  my $self = shift;

  return scalar $self->ids();
}

=head2 total_mass

Returns the sum of the mass of every map in the mapset

=cut
sub total_mass {
  my $self = shift;
  my $rvalue = 0;
  $rvalue += $_->size for $self->maps;
  $rvalue;
}

=head2 write

Write mapset information to a file.  If the method argument is a string,
then the file is opened and the maps are written to the file.  If the method
argument is a file handle, then the output is written to the file handle.

  $mapset->write('filename.maps');
  $mapset->write($fh);

=cut

sub write {
  my $self = shift;
  my $file = shift;

  my $fh = undef;
  if (ref($file) eq '') {
    open($fh, ">$file") || die("Could not open '$file' to write MapSet");
  } else {
    $fh = $file;
  }

  foreach ($self->maps) {
    $_->write($fh);
  }

  close($fh) if (ref($file) eq '');
}

=head2 read

Method read() reads data from a file name and adds those maps to the mapset.
An exception is thrown if the file cannot be opened for reading.

  $mapset->read('filename.maps');
  $mapset->read($fh);

=cut

sub read {
  my $self = shift;
  my $file = shift;

  my $fh;
  if ( ref($file) eq '') {
    open($fh, $file) || die("Could not open '$file' to read MapSet");
  } else {
    $fh = $file;
  }

  while (my $map = $self->_read_map($fh)) {
    $self->add($map);
  }

  close($fh) if ref($file) eq '';
}

#
# Private function that allows us to inherit from this, yet load different
# kinds of maps.
#

sub _read_map {
  my ($self, $fh) = @_;
  Maps->read($fh); 
}

1;
