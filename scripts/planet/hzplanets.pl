#!/usr/bin/perl -w

use strict;
use POSIX;

#idrm fields
my @fields = (0 .. 168);
my @rundes=("nro7m+10");
my $rl=10000;
my $nr=100;

#my $logamin=log10(0.72);
my $logamin=log10(0.3);
#my $logamax=log10(2.0);
my $logamax=log10(30.0);
#my $logmmin=log10(0.1);
my $logmmin=log10(10);
#my $logmmax=log10(10);
my $logmmax=log10(10);

my $inc; my $p; my $a;
my $pfile;

my $pi=4.0*atan(1);

my $mass;
my $nl=$rl;
my $nf=$nr;

for (my $r=0; $r<@rundes; $r++)
{
   my $dir = "${rundes[$r]}";
   if(! -d $dir)
   {
     mkdir($dir);
   }

#   for (my $f=0;$f<@fields;$f++)
#   {
#     my $base = "$dir/${rundes[$r]}.planets.${fields[$f]}.";
     my $base = "$dir/${rundes[$r]}.planets.";

     print "Writing $base files with $nl lines\n";

     for(my $i=0;$i<$nf;$i++)
     {
        $pfile = "$base$i";
        open(OUT,">$pfile") || die "Could not open output file $pfile\n";

        for(my $j=0;$j<$nl;$j++)
        {
  	  $a = 10**($logamin + ($logamax-$logamin)*rand());
	  $mass = 3.00374072e-6*10**($logmmin + rand()*($logmmax-$logmmin));
	  my $rnd = rand();
	  $inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/$pi;
	  #$inc = -90.0 + 180.0*rand();
	  $p = 360.0*rand();
	  print OUT "$mass $a $inc $p\n";
	} #end for nlines

	close(OUT);
      } #end for nfiles
#    } #end for fields
  } #end for rundes
