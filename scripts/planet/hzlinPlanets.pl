#!/usr/bin/perl -w

use strict;
use POSIX;

my @rundes=("hzlin");
my $rl=333;
my $nr=1000;

my $inc; my $p; my $a; my $m;
my $pfile;

my $pi=4.0*atan(1);

my $amin=0.72;
my $amax=2.0;
my $logmmin=-1;
my $logmmax=1;
my $nl=$rl;
my $nf=$nr;

for (my $r=0; $r<@rundes; $r++)
{
    my $dir = "${rundes[$r]}";

    if(! -d $dir)
    {
        mkdir($dir);
    }
  
    for (my $f=0;$f<169;$f++)
    {

        my $base = "$dir/${rundes[$r]}.planets.$f.";

        print "Writing $base files with $nl lines\n";

        for(my $i=0;$i<$nf;$i++)
        {
	    $pfile = "$base$i";
    	    open(OUT,">$pfile") || die "Could not open output file $pfile\n";

	    for(my $j=0;$j<$nl;$j++)
	    {
	        $a = $amin + rand()*($amax-$amin);
	        $m = 3.00374072e-6*10**($logmmin + rand()*($logmmax-$logmmin));
	        my $rnd = rand();
	        $inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/$pi;
	        #$inc = -90.0 + 180.0*rand();
	        $p = 360.0*rand();
	        print OUT "$m $a $inc $p\n";
	    } #end for nlines

	    close(OUT);
        } #end for nfiles
    } #end for fields
} #end for rundes
