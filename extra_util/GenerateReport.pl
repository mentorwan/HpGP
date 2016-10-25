#!/usr/bin/perl -w

## The purpose of this script is to generate report from all summary files.
## input: summary files for current folders
## output: four files. (potpos.txt, potneg.txt, realpos.txt, realneg.txt)

use strict;
use warnings;



my @array;


@array=("4023.summary.txt","4038.summary.txt","4056.summary.txt","4060.summary.txt",
        "4075.summary.txt","4080.summary.txt","4094.summary.txt","4100.summary.txt",
        "4124.summary.txt","4134.summary.txt","4135.summary.txt","4207.summary.txt",
        "4216.summary.txt","4223.summary.txt","4245.summary.txt","4270.summary.txt",
        "4281.summary.txt","4331.summary.txt","4356.summary.txt","4373.summary.txt",
        "4386.summary.txt","4396.summary.txt","4410.summary.txt","4430.summary.txt",
        "4444.summary.txt","4447.summary.txt");

my $line;
my %potpos;
my %potneg;
my %realpos;
my %realneg;
        
foreach my $currentinfile (@array)
{
	
(my $samplename= $currentinfile)=~ s/\.summary\.txt//;

print "Sample Name is ".$samplename."\n";

open(FILE, "< $currentinfile") or die "Can't open '$currentinfile':$!\n";

while ($line = <FILE>)
{
	 
   chomp $line;
   my @splitline = split /\s+/, $line;
  # print "0:".$splitline[0]."\t1:".$splitline[1]."\n";
   $potpos{$splitline[0]}{$samplename}=$splitline[1];
   $potneg{$splitline[0]}{$samplename}=$splitline[2];
   $realpos{$splitline[0]}{$samplename}=$splitline[3];
   $realneg{$splitline[0]}{$samplename}=$splitline[4];

}

#print "Value:".$potpos{$samplename}{"HP0547"}."\n";
}

my $outfile1 = "potpos.txt";
my $outfile2 = "potneg.txt";
my $outfile3 = "realpos.txt";
my $outfile4 = "realneg.txt";

writeoutput($outfile1,%potpos);
writeoutput($outfile2,%potneg);
writeoutput($outfile3,%realpos);
writeoutput($outfile4,%realneg);



sub writeoutput{

    my($outfile,%potpos)=@_;

    open(OUT, "> $outfile") or die "Can't open '$outfile':$!\n";
    
    print OUT "Gene\t";
    
    foreach my $name (sort keys %{$potpos{"HP0547"}})
    {	
	      print OUT $name."\t";
    }
    print OUT "\n";
    
    foreach my $name (sort keys %potpos)
    {
	      print OUT $name."\t";
	      foreach my $key (sort keys %{$potpos{$name}} )
	      {
		       print OUT $potpos{$name}{$key}."\t";
	      }
	  print OUT "\n";
   }

}

