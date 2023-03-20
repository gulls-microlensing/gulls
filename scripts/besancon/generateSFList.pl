#!/usr/bin/perl -w

use strict;
use POSIX;

if(@ARGV!=2)
{
    die "Usage: ./generateSFList.pl <output file> <band>\n";
}

open(OUT,">${ARGV[0]}") || die "Could not open output file";

#opendir(DIR,getcwd);

#my @files = readdir(DIR);

my @levels = ("faint","moder","bright");
my @ln = (0,1,2);
my $band=$ARGV[1];

my $nfields=27;
my @fields= (31,32,33,44,45,46,57,58,59,70,71,72,83,84,85,96,97,98,109,110,111,122,123,124,135,136,137);

for(my $j=0;$j<@levels;$j++)
{

    my $indir = "EUCLID-${band}-${levels[$j]}";

    my @files = glob "$indir/out-????.sf";
    my @f;

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
	
		print OUT "${f[0]} ${ln[$j]} $sa ${files[$i]}\n";
	    }
	}
    }
}

