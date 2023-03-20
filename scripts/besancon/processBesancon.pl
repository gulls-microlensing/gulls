#!/usr/bin/perl -w

use strict;
use POSIX;

if(@ARGV!=3 && @ARGV!=4)
{
    die "Usage: ./processBesancon.pl <input> <output extension> <input band> {<trim?>}\n";
}

open(IN,"<${ARGV[0]}") || die "Could not open input file";
open(OUT,">${ARGV[0]}${ARGV[1]}") || die "Could not open output file";

my $band=$ARGV[2];

my $trim=0;

if(@ARGV==4)
{
    $trim=1;
}

my $line;
my @data;
my @n;

@n = split('-',$ARGV[0]);

my $lo = floor(($n[-1]-1)/13)*0.25 - 0.4 - 0.125;
my $bo = floor(($n[-1]-1)%13)*0.25 - 3.2 - 0.125;

srand();

while($line = <IN>)
{
    chomp($line);
    if($line =~ m/[a-df-zA-DF-Z\*]/g) #a comment line - don't process
    {
	print OUT "$line\n";
    }
    else #a data line
    {
	@data = split(' ',$line);
	if(@data==30)
	{
	    my $r; my $riz; my $i; my $y; my $j; my $h;
	    if($band=~/^i$/i)
	    {
		$r = $data[1] + $data[0] - $data[2];
		$riz = 0.5*($data[0] + $data[1] + $data[0] - $data[2]);
		$i = $data[0];
		$y = 0.5*($data[0] + $data[3] + $data[0] - $data[2]);
		$j = $data[3]+$data[0]-$data[2];
		$h = $data[0]-$data[2];
	    }
	    elsif($band=~/^j$/i)
	    {
		$r = $data[1] + $data[0] - $data[3];
		$riz = 0.5*($data[1]+$data[2]) + $data[0] - $data[3];
		$i = $data[2] + $data[0] - $data[3];
		$y = 0.5*($data[0] + $data[2] + $data[0] - $data[3]);
		$j = $data[0];
		$h = $data[0] - $data[3];
	    }
	    elsif($band=~/^h$/i)
	    {
		$r=$data[1]+$data[0];
		$riz=0.5*($data[1]+$data[2])+$data[0]; 
		$i=$data[2]+$data[0]; 
		$y=0.5*($data[2]+$data[3])+$data[0];
		$j=$data[3]+$data[0]; 
		$h=$data[0];
	    }

	    my $ll = $lo + 0.25*rand();
	    my $bb = $bo + 0.25*rand();

	    $data[20] = $ll;
	    $data[21] = $bb;

	    print OUT "$riz $i $y $j $h ";
	    if(!$trim)
	    {
		for(my $ii=4;$ii<@data;$ii++)
		{
		    print OUT "${data[$ii]} ";
		}
	    }
	    print OUT "\n";
	}
	else
	{
	    print OUT "\n";
	}
    }
}
