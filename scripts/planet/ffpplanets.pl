#!/usr/bin/perl -w

use strict;
use POSIX;

my $amin=-1;
my $amax=2;

my $mmin=-1.5;
my $mmax=4;

my @fields = (0 .. 168);

my $nl=1000;
my $nf=40;

my $inc; my $p; my $a;
my $pfile;

#Mass function parameters
my $mpiv=95;
my $thresh=2;
my $slope=-0.73;
my $pivnorm=0.24;

my $pi=4.0*atan(1);

sub mf
{
    return $pivnorm * (10**$_[0]/$mpiv)**$slope;
}

#Important quanities
my $mmid=log10($mpiv*($thresh/$pivnorm)**(1.0/$slope));
my $da=$amax-$amin;
my $flatnorm = $thresh*($mmid-$mmin);
my $slopenorm = $pivnorm/$slope/($mpiv**$slope)/log(10) * (10**($mmax*$slope)-10**($mmid*$slope));
my $norm = $da * ( $flatnorm + $slopenorm );
my $ca=log(10)*$slope;
my $pnorm = $ca/(exp($ca*$mmax)-exp($ca*$mmid));

my $below=$flatnorm/($flatnorm+$slopenorm);

sub mofx
{
    my $x=$_[0];

    return 10**(log($ca*$x/$pnorm + exp($ca*$mmid))/$ca);
}
print "$mmid\n";
print "$flatnorm $slopenorm $norm $pnorm\n";
my $l0=&mofx(0);
my $l1=&mofx(1);

#print "$l0 $l1\n";
#exit;

my $dir="uniform";
if(! -d $dir)
{
    mkdir($dir);
}
my $rnd;
my $mass;

for (my $f=0;$f<@fields;$f++)
{

    my $base = "$dir/uniform.planets.${fields[$f]}.";

    print "Writing $base files with $nl lines\n";

    for(my $i=0;$i<$nf;$i++)
    {
	$pfile = "$base$i";
	open(OUT,">$pfile") || die "Could not open output file $pfile\n";

	for(my $j=0;$j<$nl;$j++)
	{
	    $a = 10**($amin + ($amax-$amin)*rand());
	    $rnd = rand();
	    $inc = 180*($rnd<0.5?acos(2*$rnd):-acos(2-2*$rnd))/$pi;
	    $p = 360.0*rand();

	    #mass: below or above cutoff?
	    $rnd = rand();
	    $mass = 3.00374072e-6*10**($mmin + $rnd*($mmax-$mmin));

	    print OUT "$mass $a $inc $p\n";
	} #end for nlines

	close(OUT);
    } #end for nfiles
} #end for fields
