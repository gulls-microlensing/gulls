#!/usr/bin/perl -w

use strict;
use POSIX;

my $mmin=1.0/3.0;
my $mmax=1.0/1000.0;

#idrm fields
#my @rundes=("kgrizuniform");
my @rundes=("kgrizinv");
#my @rl=(4000,4000,4000,2000,2000,2000,2000,2000,2000,2000);
my @rl=(1000);
#my @nr=(800,600,400,300,200,100,50,25,10,10);
my @nr=(20);

my $amin=log10(0.3);
my $amax=$amin+2;

#my $amin=1;
#my $amax=2;

my $inc; my $p; my $a; my $mass;
my $pfile;

my $pi=4.0*atan(1);

my $nl=$rl[0];
my $nf=$nr[0];

for (my $r=0; $r<@rundes; $r++)
{
    my $dir = "${rundes[$r]}";
    
    if(! -d $dir)
    {
        mkdir($dir);
    }

    open(FLD,"</home/penny/dmabuls/kgriz/KgrizVRIw.sources") || die "Could not open sources file\n";
    while(my $line=<FLD>)
    {
	chomp($line);
	my @data = split(' ',$line);
	if(@data!=7){next;}
	my $f = $data[0];

	my $base = "$dir/${rundes[$r]}.planets.$f.";
	
	print "Writing $base files with $nl lines\n";

	for(my $i=0;$i<$nf;$i++)
	{
	    $pfile = "$base$i";
	    if( -e $pfile){next;}
	    open(OUT,">$pfile") || die "Could not open output file $pfile\n";
	    
	    for(my $j=0;$j<$nl;$j++)
	    {
		$a = 10**($amin + ($amax-$amin)*rand());
		$mass = 3.00374072e-6/($mmin + rand()*($mmax-$mmin));
		#$a = 1+int(3*rand())/2;
		my $rnd = rand();
		$inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/$pi;
		#$inc = -90.0 + 180.0*rand();
		$p = 360.0*rand();
		print OUT "$mass $a $inc $p\n";
	    } #end for nlines

	    close(OUT);
	} #end for nfiles
    } #end for fields
    close(FLD);
} #end for rundes

