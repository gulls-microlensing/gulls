#!/usr/bin/perl -w

use strict;
use POSIX;

#massgrid 

my $name="hzgrid";
my $nsub=500;
my $nrep=10;
my $mearth=3.00374072e-6;
my $pi=4.0*atan(1);
my $logmp=0;

my $mmin=log10(0.4);
my $mmax=log10(1.8);
my $nm=14;

my $amin=-0.5;
my $amax=0.5;
my $na=10;

if(! -d "$name/")
{
    mkdir("$name");
}
chdir("$name/");

for(my $subrun=0;$subrun<$nsub;$subrun++)
{
    print "run $subrun\n";
    for(my $field=0;$field<169;$field++)
    {
	my $fname = "$name.planets.$field.$subrun";
	open(FILE,">$fname") || die "Could not open file $fname";
	for(my $n=0;$n<$nrep;$n++)
	{
	    for(my $m=0;$m<$nm;$m++)
	    {
		#log m = -2 to 1 every 0.25
		#13
		my $mass = int(100*$logmp) + 10**($mmin+$m*($mmax-$mmin)/($nm-1.0));
		for(my $aa=0;$aa<$na;$aa++)
		{
		    #17
		    #13*17=221
		    #log a = -1 to 1 every 0.125
		    $a = 10**($amin+$aa*($amax-$amin)/($na-1.0));
		    my $rnd = rand();
		    my $inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/$pi;
		    my $phase = 360*rand();
		    print FILE "$mass $a $inc $phase\n";
		}
	    }
	}
	close(FILE);
    }
}
