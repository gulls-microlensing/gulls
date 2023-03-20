#!/usr/bin/perl -w

use strict;

my @mp=("-10","+00","+10","+20","+30","+40");
my @band=("H", "J", "Y", "V"); #bands
my @princobs=(0,1,2,3);
my @rl=(3000,2000,1000,1000,1000,1000); #lenghth of runs
my @nr=(3000,2000,1000,1000,1000,1000); #number of runs

my $amin=-1.52;
my $amax=3+$amin;
my $da=$amax-$amin;
my $a; my $inc; my $p;
my $pfile;

for (my $m=0; $m<@mp; $m++)
{
  my $mm=int($mp[$m])/10.0;
  my $mass = 3.00374072e-6*10**($mm); #mp is a factor of ten larger
  my $nl=$rl[$m];
  my $nf=$nr[$m];

  for (my $b=0; $b<@band; $b++)
  {
    my $dir = "b${band[$b]}m${mp[$m]}";

    if(! -d $dir)
    {
      mkdir($dir);
    }

    my $base = "$dir/b${band[$b]}m${mp[$m]}planets.";

    print "Writing $base files with $nl lines\n";

    for(my $i=0;$i<$nf;$i++)
    {
	$pfile = "$base$i";
	open(OUT,">$pfile") || die "Could not open output file $pfile\n";

	for(my $j=0;$j<$nl;$j++)
	{
	    $a = 10**($amin+$da*rand());
	    my $rnd = rand();
	    $inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/pi;
	    $p = 360.0*rand();
	    print OUT "$mass $a $inc $p\n";
	} #end for nlines

	close(OUT);
    } #end for nfiles
  } #end for b
} #end for m
