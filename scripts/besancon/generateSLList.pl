#!/usr/bin/perl -w

use strict;
use POSIX;

if(@ARGV!=2)
{
    die "Usage: ./generateSLList.pl <output file> <input directory>\n";
}

open(OUT,">${ARGV[0]}") || die "Could not open output file";

#opendir(DIR,getcwd);

#my @files = readdir(DIR);

my $indir=$ARGV[1];

$indir=~s£/$££;

my @files = glob "$indir/out-????.sl";
my @f;

my $nfields=27;
my @fields= (31,32,33,44,45,46,57,58,59,70,71,72,83,84,85,96,97,98,109,110,111,122,123,124,135,136,137);

#grab the area
open(IN,$files[0]) || die "Could not open first file\n";

my $line;
my $sa=0;
my @data;

while($line=<IN>)
{
    chomp($line);
    if($line =~ /Solid angle/)
    {
	@data = split(' ',$line);
	my $target=-1;
	for(my $i=0;$i<@data;$i++)
	{
	    if($data[$i] =~ /angle/)
	    {
		$target=$i+1;
	    }
	    if($i==$target)
	    {
		$sa = $data[$target];
		last;
	    }
	}
    }
    if($sa>0) 
    {
	last;
    }
}

for(my $i=0;$i<@files;$i++)
{

    for(my $jj=0;$jj<$nfields;$jj++)
    {
	if($i==$fields[$jj])
	{
	    my @n = split('-',$files[$i]);
	    @f = split('\.',$n[-1]);
	    
	    $f[0]-=1;
	    
	    my $l = floor($f[0]/13)*0.25 - 0.4;
	    my $b = floor($f[0]%13)*0.25 - 3.2;

	    print OUT "${f[0]} $l $b $sa $files[$i]\n";
	}
    }
}

