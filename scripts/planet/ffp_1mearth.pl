#!/usr/bin/perl -w

use strict;
use POSIX;

#massgrid 

my $name="ffp_1mearth";
my $nsub=50;
my $nrep=10000;
my $mearth=3.00374072e-6;
my $pi=4.0*atan(1);

if(! -d "$name/")
{
    mkdir("$name");
}
chdir("$name/");

for(my $subrun=0;$subrun<$nsub;$subrun++)
{
    for(my $field=0;$field<169;$field++)
    {
	my $fname = "$name.planets.$field.$subrun";
	open(FILE,">$fname") || die "Could not open file $fname";
	for(my $n=0;$n<$nrep;$n++)
	{
	    #log m = -2 to 1 every 0.25
	    #13
	    my $mass = $mearth*10**1;
	    $a = 1;
	    my $rnd = rand();
	    my $inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/$pi;
	    my $phase = 360*rand();
	    print FILE "$mass $a $inc $phase\n";
	}
	close(FILE);
    }
}
