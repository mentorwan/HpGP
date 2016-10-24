#!/usr/bin/perl -w

## The purpose of this script is to generate location for each motification motif
## Input: motif.gff file from PacB output and methylation motif list
## Output: a table to get circos plot


use strict;
use warnings;
use Bio::SeqIO;
use File::Basename;

@ARGV == 4 or die 
"usage: $0 motif.gff motif_list.txt polished.fasta gene.fasta\n";

#system("module load blast");

my ($infile, $motiffile, $fasta1, $fasta2) = @ARGV;


##Reading motif pattern from motif file. If one motif has multiple motify types (m6A,m4C), it will report an issue.

open(FILE1, "< $motiffile") or die "Can't open '$motiffile':$!\n";

my $line1;
my %motif;
my $count=0;
while ($line1 = <FILE1>)
{
	chomp $line1;	
	$line1 =~ s/\s+$//;   #remove whitespace after string.
	$motif{$line1} = 1;
	$count++;	
}

#remove suffix and get directory
(my $directory= basename($infile))=~ s/\.[^.]+$//;;

#print "Total motif count is ".$count."\n";

my %sequences;
my $seqio = Bio::SeqIO->new(-file => $fasta2);
while(my $seqobj = $seqio->next_seq )
{
	 my $id = $seqobj->display_id;
	 my $seq = $seqobj->seq;
	 $sequences{$id}=$seq;
	 print $id."\n";
	 my $seqout=Bio::SeqIO->new(-file => ">$directory/temp.fasta", -format => "fasta");
	 $seqout->write_seq($seqobj);



#run the blast first.
              
my $result = `blastn -query $fasta1 -subject $directory/temp.fasta -outfmt 6 -max_target_seqs 1`;

print "##".$result;

if ( $result eq '')
{
	print "No alignment found.\n";
	next;
}

my @resultarray = split /\t/, $result;

#print $resultarray[0]."\t".$resultarray[6]."\t".$resultarray[7]."\n";

##if ( $resultarray[0] !~ /unitig_0/ )
##{
##	print "BLAST is failed.\n";
##	exit;
##}

if ($resultarray[6] > $resultarray[7] )
{
	print "The alignment is reversed. Need to check.\n";
	exit;
}



##Read input file from GFF file

open(FILE, "< $infile") or die "Can't open '$infile':$!\n";

my $line ;

#print "Contig\tStart\tEnd\tStrand\tMotif\tType\tQVscore\n";
while ($line = <FILE>){

  chomp $line;
  my @splitarray;
  my $substring;
  my $QVscore;
  my $contignumber;
  
  
  if ($line !~ /motif=/) { next;}
  else {
  	my @splitarray = split /\t/, $line;
  
       	
  	if ($splitarray[0] eq $resultarray[0] and $splitarray[2] ne '.' and  $splitarray[2] ne 'modified_base' and $splitarray[3] >= $resultarray[6] and $splitarray[3] <= $resultarray[7] )  ## only select largest contig, supposely the first contig
  	{ 
  		   
  		  #print $splitarray[3]."\t". $resultarray[6]."\t".$resultarray[7]."\n";
  		  ($substring) = $splitarray[8]=~/motif=(.*);coverage/;  ##find motif
  		  
  		  ($QVscore) = $splitarray[8]=~/identificationQv=(.*)$/;  #find QVscore
  		  
  		  ($contignumber) = $splitarray[0] =~ /unitig_(\d+)|quiver/;
  		  
  		  if( exists $motif{$substring} )
  	    {
         	print "0\t$splitarray[3]\t$splitarray[4]\t$splitarray[6]\t$substring\t$splitarray[2]\t$QVscore\n";
  	    }
  	}
  }
}

close FILE;
}