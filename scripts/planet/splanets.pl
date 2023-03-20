#!/usr/bin/perl -w

use strict;
use POSIX;

#idrm fields
my @fields = (0 .. 168);
my @mp=("-250","-225","-200","-175","-150","-125","-100");
#my @mp=("+00", "+20");
my @rundes=("l");
#my @rl=(4000,4000,4000,2000,2000,2000,2000,2000,2000,2000);
#my @rl=(3333, 3333, 3333, 3333);
my @rl=(1000, 1000, 1000, 1000, 1000, 1000, 1000);
#my @nr=(800,600,400,300,200,100,50,25,10,10);
my @nr=(15, 15, 15, 15, 15, 15, 15);

#my $amin=log10(2);
#my $amax=$amin+0;

my $amin=-0.25;
my $astep=0.125;
my $amax=2;

my $inc; my $p;
my $pfile;

my $pi=4.0*atan(1);

for (my $m=0; $m<@mp; $m++)
{

  my $mm=int($mp[$m])/100.0;
  my $mass = 3.00374072e-6*10**($mm); #mp is a factor of ten larger
  my $nl=$rl[$m];
  my $nf=$nr[$m];

  for(my $a=0; $a<($amax-$amin)/$astep; $a++)
  {
    my $aa = 10**($amin+$a*$astep);
    for (my $r=0; $r<@rundes; $r++)
    {
      my $dir = "${rundes[$r]}m${mp[$m]}a$a";

      if(! -d $dir)
      {
        mkdir($dir);
      }

      for (my $f=0;$f<@fields;$f++)
      {

	  my $base = "$dir/${rundes[$r]}m${mp[$m]}a$a.planets.${fields[$f]}.";

	  print "Writing $base files with $nl lines\n";

	  for(my $i=0;$i<$nf;$i++)
	  {
	      $pfile = "$base$i";
	      open(OUT,">$pfile") || die "Could not open output file $pfile\n";

	      for(my $j=0;$j<$nl;$j++)
	      {
#		  $a = 10**($amin + ($amax-$amin)*rand());
		  my $rnd = rand();
		  $inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/$pi;
		  #$inc = -90.0 + 180.0*rand();
		  $p = 360.0*rand();
		  print OUT "$mass $aa $inc $p\n";
	      } #end for nlines

	      close(OUT);
	  } #end for nfiles
      } #end for fields
    } #end for rundes
  } #end for a
} #end for m
