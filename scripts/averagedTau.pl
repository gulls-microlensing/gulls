#!/usr/bin/perl -w

use strict;
use POSIX;

if(@ARGV!=2&&@ARGV!=3)
{
    die "Usage: ./calcRates.pl <Source file> <Lens file> {<nfilters>=5}\n";
}

my $sname = $ARGV[0];
my $lname = $ARGV[1];

#Less commonly changed user options
my $nfilters = 5;
if(@ARGV==3) {$nfilters=$ARGV[2];}
my $mulcol = 0 + $nfilters;
my $mubcol = 1 + $nfilters;
my $distcol = 20 + $nfilters;
my $masscol = 12 + $nfilters;
my $carea = 0.25**2; #Area of a line of sight

my $uucol = 4 + $nfilters;
my $vvcol = 5 + $nfilters;
my $wwcol = 6 + $nfilters;
my $lcol = 17 + $nfilters;
my $bcol = 18 + $nfilters;

#constants
my $rEsun = 2.85412; 
#multiply by this factor to convert from mas yr-1 to km s-1 / kpc
my $mAUYRSEC = 1.495978707e8 / (365.25 * 86400.0);
my $pi = 3.1415926535897932384626433832795;

#Variables
my $sarea=0.0002; #Area used to generate source file
my $larea=0.0008; #Area used to generate lens area
#my $larea=0.00020; #Area used to generate lens area

my @lmass; my @ldist; my @lmul; my @lmub;
my $tausum; my $tEsum; my $wsum; my $nod;

my @efftE = (-0.9, -0.7, -0.5, -0.3, -0.1, 0.1, 0.3, 0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9, 2.1, 2.3);
#my $nsub=78; my @eff = (0, 0, 0.003273, 0.007095, 0.011073, 0.016599, 0.022921, 0.030432, 0.036723, 0.04443, 0.052923, 0.056168, 0.061749, 0.056394, 0.048691, 0.030679, 0.004935); #MOA field 4
#my $nsub=65; my @eff = (0, 0, 0.004388, 0.010803, 0.014527, 0.020526, 0.026112, 0.034137, 0.039902, 0.04291, 0.051266, 0.057265, 0.060719, 0.06245, 0.050169, 0.033799, 0.006419); #MOA field 5
my $nsub=11; my @eff = (0, 0, 0.002482, 0.006764, 0.010261, 0.015368, 0.021227, 0.028904, 0.031804, 0.039, 0.047719, 0.0531, 0.058713, 0.052162, 0.039349, 0.019411, 0.002302); #MOA field 6
#my $nsub=78; my @eff = (0, 0, 0.003961, 0.009205, 0.012539, 0.018845, 0.026169, 0.035786, 0.045184, 0.053996, 0.062324, 0.070292, 0.075805, 0.072012, 0.061701, 0.036734, 0.00474); #MOA field 8
#my $nsub=79; my @eff = (0, 0, 0.005647, 0.01094, 0.017164, 0.02416, 0.02975, 0.039277, 0.048125, 0.055136, 0.062355, 0.068945, 0.071557, 0.069448, 0.062441, 0.040306, 0.007455); #MOA field 9
#my $nsub=70; my @eff = (0, 0, 0.002863, 0.006496, 0.009337, 0.014661, 0.021471, 0.027825, 0.035548, 0.042866, 0.051088, 0.057712, 0.059359, 0.05683, 0.048392, 0.025402, 0.002706); #MOA field 10
#my $nsub=79; my @eff = (0, 0, 0.003467, 0.007747, 0.011062, 0.016986, 0.023694, 0.032517, 0.041371, 0.050045, 0.060234, 0.064146, 0.070571, 0.067371, 0.058694, 0.029789, 0.004014); #MOA field 14


#print a header line
print "#field <tE/day> <tau> <nev/yr-1>\n";

print "$sname $lname\n";

#Preload the lens data - assume it is sorted by distance
open(LENS,"<$lname") || die "Could not open lens file $lname\n";

#Get the 
    
my $lline;
my @ldata;

#Clear the lens arrays
@lmul = ();
@lmub = ();
@ldist = ();
@lmass = ();

#Clear the accumulators
$tausum = 0;
$tEsum = 0;
$wsum = 0;
$nod=0;

my $timescalesum=0;

while($lline = <LENS>)
{
    chomp($lline);
    
    #Get the lens file area
    if($lline =~ /Solid angle/)
    {
	@ldata = split(' ',$lline);
	for(my $i=0;$i<@ldata;$i++)
	{
	    #print "${ldata[$i]}\n";
	    if($ldata[$i] =~ /angle/)
	    {
		$larea = $ldata[$i+1];
		last;
	    }
	}
    }
    
    #Get the data of all lenses
    if(!($lline =~ /[a-df-zA-DF-Z*]/))
    {
	@ldata = split(' ',$lline);
	
	if(@ldata==26+$nfilters)
	{
	    push(@ldist, $ldata[$distcol]);
	    push(@lmass, $ldata[$masscol]);
#	    if($nfilters==4)
#	    {
		my $dl=$ldata[$distcol];
		my $ul=0.2109454*($ldata[$uucol]-10.3)/$dl; 
		my $vl=0.2109454*($ldata[$vvcol]-232.8)/$dl; 
		my $wl=0.2109454*($ldata[$wwcol]-5.9)/$dl;
		my $ll=$ldata[$lcol]*$pi/180.0; 
		my $bl=$ldata[$bcol]*$pi/180.0;

#		push(@lmul,$ul*sin($ll)*cos($bl) + $vl*cos($ll)*cos($bl) + $wl*sin($bl));
#		push(@lmub,$ul*sin($bl) + $vl*sin($ll)*cos($bl) + $wl*cos($bl));
		
#	    }
#	    else
#	    {
		push(@lmul, $ldata[$mulcol]);
		push(@lmub, $ldata[$mubcol]);
#	    }

	    #print "${lmul[-1]} ${lmub[-1]} ${ldist[-1]} ${lmass[-1]}\n";
	}
    }
}

close(LENS);

#Now load the sources
open(SRC,"<$sname") || die "Could not open source file $sname\n";

my $sline;
my @sdata;

my $tau_average=0;
my $tau_rw=0; #rate weighted average
my $tau_tE=0; #rate weighted average
my $nsrc=0;
my $rwsrc=0;
my $meantEsum=0;
my $rEsum=0;
my $moasum=0;
my $moaeff;

#Read in the sources and process
while($sline = <SRC>)
{
    chomp($sline);
    
    #Get the source file area
    if($sline =~ /Solid angle/)
    {
	@sdata = split(' ',$sline);
	for(my $i=0;$i<@sdata;$i++)
	{
	    if($sdata[$i] =~ /angle/)
	    {
		$sarea = $sdata[$i+1];
		last;
	    }
	}
    }
    if(!($sline =~ /[a-df-zA-DF-Z*]/))
    {
	@sdata = split(' ',$sline);
	
	if(@sdata==26+$nfilters)
	{
	    #Basic data
	    my $smul;
	    my $smub;
#	    if($nfilters==4)
#	    {
		my $ds=$sdata[$distcol];
		my $us=0.2109454*($sdata[$uucol]-10.3)/$ds;
		my $vs=0.2109454*($sdata[$vvcol]-232.8)/$ds; 
		my $ws=0.2109454*($sdata[$wwcol]-5.9)/$ds;
		my $ls=$sdata[$lcol]*$pi/180.0; 
		my $bs=$sdata[$bcol]*$pi/180.0;

#		$smul = $us*sin($ls)*cos($bs) + $vs*cos($ls)*cos($bs) + $ws*sin($bs);
#		$smub = $us*sin($bs) + $vs*sin($ls)*cos($bs) + $ws*cos($bs);
		



#	    }
#	    else
#	    {
		$smul = $sdata[$mulcol];
		$smub = $sdata[$mubcol];
#	    }
	    my $sdist = $sdata[$distcol];
	    
	    my $tau = 0; #reset the optical depth accumulator
	    my $rwtau = 0; #reset the optical depth accumulator
	    my $rw=0;
	    my $nlens=0;
	    my $tEmean=0;

	    #Now calculate things for each source lens pair
	    
	    for(my $l=0; $l<@ldist; $l++)
	    {
		if($sdist<=$ldist[$l])
		{
		    #Lens is further than the source
		    last;
		}
		
		my $x = $ldist[$l]/$sdist;
		my $rE = $rEsun * sqrt($lmass[$l] * $sdist * (1-$x) * $x);
		my $thE = $rE / $ldist[$l];
		my $murel = sqrt(($smul-$lmul[$l])**2 + ($smub-$lmub[$l])**2) + 1.0e-50;
		my $vt = $murel * $ldist[$l] * $mAUYRSEC;
		my $tE = 365.25 * $thE / $murel;
		my $w = 2 * $thE * $murel;

		my $logtE=log10($tE);

		if($logtE>=2||$logtE<-1){$moaeff=0;}
		else
		{
		    my $idx=int(20*($logtE+1)/4.0);
		    #my $mindiff=1e30;
		    #my $diff;
		    #for(my $tt=0;$tt<@efftE;$tt++)
		    #{
	#		$diff = abs($logtE-$efftE[$tt]);
	#		if($diff<$mindiff)
	#		{
	#		    $mindiff=$diff;
	#		    $idx=$tt;
	#		}
	#	    }
		    $moaeff = $eff[$idx];
		}

		$tau += $thE**2;
		$rwtau += $thE**2;
#		$tEtau += $tE * $thE**2;
		$tEmean += $tE;
		$timescalesum += $tE*$w;
		$rEsum += $rE*$w;
		
		#accumulate the averages and weightings
		$wsum += $w;
		$moasum += $moaeff * $w;
		$rw += $w;
		$nlens++;
	    } #end for each lens

	    $tau *= $pi / (1000**2 * $larea*3600**2);
	    $rwtau *= $pi / (1000**2 * $larea*3600**2);
	    if($nlens>0) {$tEmean /= $nlens;}

	    if($sdist==8) {print "tau(8kpc) = $tau\n";}

	    $tau_average+=$tau;
	    $tau_rw+=$rw*$tau; #rate weighted average
	    if($nlens>0)
	    {
		$tau_tE+=$tau*$rw*$tEmean;
		$meantEsum+=$rw*$tEmean;
	    }
	    $nsrc++;

	    #print "tau: $sdist $tau $rwtau\n";

	} #end if sdata valid
    } #end if ! header line
} #end for each source

close(SRC);

$tau_average /= $nsrc;
$tau_rw /= $wsum;
$tau_tE /= $meantEsum;
my $meantE = $timescalesum/$wsum;
my $meanrE = $rEsum/$wsum;

print "Field averaged: tau_average\ttau_rw\ttau_rw/tE\n";
print "Field averaged: $tau_average\t$tau_rw\t$tau_tE\n";
print "Mean timescale: $meantE\n";
print "Mean rE: $meanrE\n";

#Now work out the final values

#and the event rate
my $nev = $wsum * 1e-6 / ($larea*3600*3600) * (0.25 * 0.25 / $sarea);
my $nev_per_src = $nev / (0.25*0.25*$nsrc/$sarea);
my $nev_per_deg2 = $nev/(0.25*0.25);

my $nevmoa = $moasum * 1e-6 / ($larea*3600*3600) * (8096*0.58/3600 * 10240*0.58/3600 / $sarea) * (596.0/365.25) * ($nsub/80.0);
# is the sum of individual swept areas / the area of the lens file, 
#multiplied by the ratio of cell to source file areas

#Print the output
print "rate: Rate_per_tile Rate_per_deg^2 Rate_per_src\n";
print "rate: $nev $nev_per_deg2 $nev_per_src\n";
print "moarate: $nevmoa\n";

