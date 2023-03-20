#!/usr/bin/perl -w

use strict;

if (@ARGV != 2 && @ARGV != 3)
{
    print STDERR "Usage: ./applycovfac.sh <input> <covfac> {<ratecol>}\n";
    die "\tOutput file is <input>.covfac\n";
    
}
   
my $inname=$ARGV[0];
my $cfname=$ARGV[1];
my $ratecol=-1;

if(@ARGV==3) 
{
    $ratecol=$ARGV[2];
}


#Name the output file
if($cfname !~ m/\//)
{
    if( -e $cfname)
    {
	$cfname="./$cfname";
    }
    else
    {
	$cfname=$ENV{"HOME"} . "/dmabuls/rates/$cfname";
    }
}
$cfname =~ /.*\/(.+)\.covfac$/ || die "covfac file must end in .covfac\n";
my $cfroot=$1;
my $outname="$inname.${cfroot}_covfac";

#print "$inname $cfroot $outname $cfname\n";
#exit(1);

#Open the files

open(IN,"<$inname") || die "Could not open input file ($inname)\n";
open(CF,"<$cfname") || die "Could not open covfac file ($cfname)\n";
open(OUT,">$outname") || die "Could not open output file ($outname)\n";

my $line;
my @data;

#Load the covering fractions
my @covfac = (0) x 169;
while($line=<CF>)
{
    chomp($line);
    @data = split(' ',$line);
    
    if(@data>=2)
    {
	$covfac[$data[0]] = $data[1];
    }
}

#Parse the input
while($line=<IN>)
{
    chomp($line);
    @data = split(' ',$line);

    if(@data>$ratecol)
    {
	my $field=$data[2];
	#print "$covfac[$field] $field $data[$ratecol]\n";
	if($covfac[$field]>0)
	{
	    $data[$ratecol] *= $covfac[$field];
	    print OUT "@data\n";
	}
    }
}
