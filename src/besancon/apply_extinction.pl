#!/usr/bin/perl -w

use strict;
use POSIX;

if(@ARGV!=1)
{
    die "Usage: ./apply_extinction.pl <input>\n";
}

#parse a besancon file to convert colors to magnitudes, apply extinction and rearrange the columns

my $Mbolsun = 4.7554;
my $Teffsun = 5771.8;
my $R0 = 8.5;

my $sourcelimit = 23.0;

my $extfile = "marshall_extinction.txt";

my $infile = $ARGV[0];
my $outfile = "kgriz/$infile.KgrizVRIw";
my $headfile = "$infile.header";

if(!($infile =~ /^l([0-9\.+-]+)_b([0-9\.+-]+)_(.+)$/))
{
    die "Unrecognized filename format - expect l<l>_b<b>_<type>\n";
}
my $l = $1;
my $b = $2;
my $type = $3;

print "$l $b $type\n";

my $sourcetype=0;
if($type eq "source")
{
    $sourcetype=1;
}

open(IN,"<$infile") || die "Could not open input file $infile\n";
open(EXT,"<$extfile") || die "Could not open extinction file $extfile\n";
open(OUT,">$outfile") || die "Could not open output file $outfile\n";

#Calculate cartesian position and ra,dec
my $pi = 4*atan(1);
my $cosb = cos($b*$pi/180.0);
my $sinb = sin($b*$pi/180.0);
my $cosl = cos($l*$pi/180.0);
my $sinl = sin($l*$pi/180.0);

#my $sindngp=0.455983776;
#my $cosdngp=0.889988087;
#my $sinangp=-0.222560701;
#my $cosangp=-0.974918834;
#my $rangp=3.36603292;
#my $ra0=4.936838322;
#my $l0=0.574736922;
#my $sinlml0=sin($l-$l0); 
#my $coslml0=cos($l-$l0);
#my $sind = $sinb*$sindngp + $cosb*$cosdngp*$sinlml0;
#my $cosd = sqrt(1.0-$sind*$sind);
#my $cosama0 = $coslml0*$cosb/$cosd;
#my $sinama0 = (-$sinb*$cosdngp + $cosb*$sindngp*$sinlml0)/$cosd;  
#my $ra = $ra0 + atan2($sinama0,$cosama0);
#$ra *= 180.0/$pi;
#my $dec = asin($sind);
#$dec *= 180.0/$pi;

my @jgal = ( [-0.054875539726,-0.873437108010,-0.483834985808],
	     [ 0.494109453312,-0.444829589425, 0.746982251810],
	     [-0.867666135858,-0.198076386122, 0.455983795705]);

#  Spherical to Cartesian
my $radl = $l*$pi/180.0;
my $radb = $b*$pi/180.0;
my $r_ = 1.0;
my @pos = ( cos($radl)*cos($radb), sin($radl)*cos($radb), sin($radb) );

#  Rotate to equatorial coordinates 
my @pos1 = (0 x 3);
for (my $i = 0; $i < 3; $i++) 
{
    $pos1[$i] = $pos[0]*$jgal[0][$i] + $pos[1]*$jgal[1][$i] 
	+ $pos[2]*$jgal[2][$i];
}

#  Cartesian to Spherical
my $ra = atan2($pos1[1], $pos1[0]);
if ($ra < 0) {$ra = $ra + 2*$pi;}
my $rxy2 = $pos1[0]**2 + $pos1[1]**2;
my $rxy = sqrt($rxy2);
my $dec = atan2($pos1[2],$rxy);
$ra *= 180.0/$pi;
$dec *= 180.0/$pi;

print "$ra $dec\n";

#Read in the extinction
my @dext;
my @ext;
my $found=0;
while(my $line=<EXT>)
{
    chomp($line);
    my @data = split(' ',$line);
    if(@data>=5)
    {
	my $lext = $data[0];
	my $bext = $data[1];
	if($lext==$l && $bext==$b)
	{
	    for(my $i=3;$i<@data;$i+=2)
	    {
		push(@dext,$data[$i]);
		push(@ext,$data[$i+1]);
		#print "$dext[-1] $ext[-1]\n";
	    }
	    $found++;
	}
    }
}
if($found!=1)
{
    die "The sightline was found $found times in the extinction file. Sighlines must be snapped onto the extinction grid and must only appear in the extinction file once\n";
}
close(EXT);

#Now we can process the file

my $key="     u      g      r       i       z  Kepler      V      R      I wideVR    Z087   J    H    K    W149 ";
$key="$key    mul     mub      Vr      UU      VV      WW ";
$key="$key    Mv CL  Typ  Teff  logg Age   Mass   Mbol   Radius [Fe/H] ";
$key="$key        l         b        ra       dec ";
$key="$key  Dist       X       Y       Z     A_K [a/Fe]";
my $keyadded=0;

my $lastD=0;

#my @AoAk_ugriz = (15.1555148688851, 11.6372437823174, 7.93216929311387, 
#		  5.90292118482206, 4.20796785004286); #Megacam
my @AoAk_ugriz = (16.5454766990806, 11.9846410966622, 7.97321876609085, 
		  6.00592458395235, 3.99468494156202); #SDSS
my @AoAk_ZJHKW = (3.99468494156202, 2.50672923909721, 1.57226125053936, 0.977308619864285, 1.57226125053936); #WFIRST

while(my $line=<IN>)
{
    chomp($line);
    my @data = split(' ',$line);
    
#Bad line:   5.144 -3.50  4 2.10 4.430  3.90  113.50 -0.297 -0.472 -0.370 -0.276 10.728   -0.076   -0.032  -17.59  -17.41  -18.76   -7.77  0.47   0.75000  -0.50000  0.000 -3.219
#            0.827 12.70  5 7.50 3.511  5.01  7 0.23  2.514  1.181  1.283  0.554 20.574   -1.732   -0.099    3.94    4.79  -67.81   -3.92 -0.29   0.75000  -0.50000  0.000 11.044
    if(@data!=24 || $line =~ /[a-zA-Z\*]/)
    {
	if(!($line =~ /#/))
	{
	    if(@data==23 && $data[6]>100)
	    {
		#>10 msun causes a rare problem
		$data[6] =~ /^(\d+)(\d\d\.\d\d)$/;
		my $stelpop=$1;
		my $stelmass=$2;
		#print STDERR "Bad data: @data\n";
		splice(@data,6,0,$stelpop);
		$data[7] = $stelmass;
		#print STDERR "Good data: @data\n";		
	    }
	    else
	    {
		print STDERR "$infile Bad line: $line\n";
	    }
	}
	else
	{
	    print OUT "$line\n";
	}
	next;
    }

    if(!$keyadded)
    {
	print OUT "#   l = $l   b = $b   type = $type\n";
	print OUT "#  ugriz extinction: ${AoAk_ugriz[0]} ${AoAk_ugriz[1]} ${AoAk_ugriz[2]} ${AoAk_ugriz[3]} ${AoAk_ugriz[4]}\n";
	print OUT "#$key\n";
	$keyadded=1;
    }

    my ($D, $Mv, $CL, $Typ, $Teff, $logg, $Age, $Mass, $umg, $gmr, $rmi, $imz, $i_, $mul, $mub, $Vr, $UU, $VV, $WW, $FEH, $l, $b, $A_K, $Mbol) = @data;
    my ($up,$gp,$rp,$ip,$zp); #MEGACAM magnitudes AB
    my ($us,$gs,$rs,$is,$zs); #SDSS magnitudes AB
    my ($kp,$V,$R,$I,$wVR); #Kepler, Johnson & DECam wideVR
    my ($Z087,$J,$H,$K,$W149); #WFIRST, AB system

    #Calculate the extinction - assumes sorted by distance
    #print "$D $dext[$lastD]\n";
    while($lastD<@dext && $D>$dext[$lastD])
    {
	$lastD++;
    }
    if($lastD==@dext)
    {
	$A_K = $ext[$lastD-1];
	#print "$D ${dext[$lastD-1]} inf $A_K\n";
    }
    elsif($lastD==0)
    {
	$A_K = $D/$dext[0]*$ext[0];
	#print "$D 0 ${dext[$lastD]} $A_K\n";
    }
    else
    {
	$A_K = $ext[$lastD-1] + ($D-$dext[$lastD-1])*($ext[$lastD]-$ext[$lastD-1])/($dext[$lastD]-$dext[$lastD-1]);
	#print "$D ${dext[$lastD-1]} ${dext[$lastD]} $A_K\n";
    }

    #Convert magnitudes from the Megacam AB system to the desired system

    #Megacam magnitudes
    $ip = $i_;
    $zp = $ip - $imz;
    $rp = $ip + $rmi;
    $gp = $rp + $gmr;
    $up = $gp + $umg;

    #Convert to sdss magnitudes
    $us = $up + 0.181 * ($up - $gp);
    $gs = $gp + 0.195 * ($gp - $rp);
    $rs = $rp + 0.011 * ($gp - $rp);
    #$is = $ip + 0.079 * ($rp - $ip);  #(old filter)
    $is = $ip + 0.044 * ($gp - $ip);  #(old filter)
    $zs = $zp - 0.099 * ($ip - $zp);

    #Bilir et al 2008 (assume WFIRST filters same as JHK_AB)
    $J = $gs - 1.379 * ($gs-$rs) - 1.702 * ($rs-$is) - 0.518 + 0.91;
    $H = $gs - 1.849 * ($gs-$rs) - 1.536 * ($rs-$is) - 0.666 + 1.39;
    $K = $gs - 1.907 * ($gs-$rs) - 1.654 * ($rs-$is) - 0.684 + 1.85;

    $Z087 = $zs;
    $W149 = ($J+$H+$K)/3.0;
	
    #Assume that the DECam system is the same as SDSS - may want to change later
    
    #Add extinction before calculating other magnitudes
    $us += $A_K*$AoAk_ugriz[0];
    $gs += $A_K*$AoAk_ugriz[1];
    $rs += $A_K*$AoAk_ugriz[2];
    $is += $A_K*$AoAk_ugriz[3];
    $zs += $A_K*$AoAk_ugriz[4];

    $Z087 += $A_K*$AoAk_ZJHKW[0];
    $J    += $A_K*$AoAk_ZJHKW[1];
    $H    += $A_K*$AoAk_ZJHKW[2];
    $K    += $A_K*$AoAk_ZJHKW[3];
    $W149 += $A_K*$AoAk_ZJHKW[4];

    $kp = ($gs-$rs>0.8 ? 0.1*$gs + 0.9*$rs : 0.2*$gs + 0.8*$rs);

    #Lupton (2005)
    $V = $gs - 0.5784 * ($gs - $rs) - 0.0038;  #sigma = 0.0054
    $R = $rs - 0.2936 * ($rs - $is) - 0.1439;  #sigma = 0.0072
    $I = $is - 0.3780 * ($is - $zs)  -0.3974;  #sigma = 0.0063 
    $wVR = 0.5*($V+$R);



    if($sourcetype && $is>$sourcelimit)
    {
	#Source is not bright enough
	next;
    }

    #Convert proper motions of mas/yr from "/cen
    $mul *= 10.0;
    $mub *= 10.0;

    #Calculate the radius
    my $lum = 10**(-0.4*($Mbol-$Mbolsun));
    $Teff = 10**$Teff;
    my $Radius = sqrt($lum)/($Teff/$Teffsun)**2;

    #Different for white dwarfs, use:
    if($Typ==9)
    {
	$Radius = 0.01 * $Mass**(1.0/3.0);
    }

    #Correct things that have erroneous 0 bolometric magnitude
    if($Mbol==0&&$CL==5&&$Mass<0.15)
    {
	$Radius = 10**(0.111 + 4.0/3.0*log10($Mass));
    }
   

    #Cartesian position
    my $X = $D*$cosb*$cosl;
    my $Y = $D*$cosb*$sinl;
    my $Z = $D*$sinb;
    

    my $alpha = 0.0; #[Alpha/Fe]


    #Correct U,V,W and RV,mul,mub - will be in c++ code after this
    #if($X>8.5){
	#$VV = -452-$VV;
    #}

    #my $commstr = "./kgrizcoordconvert $l $b $D $UU $VV $WW";
    #my $tmpstr = `$commstr`;
    #my @vmu = split(" ",$tmpstr);
    #$Vr = $vmu[0];
    #$mul = $vmu[3];
    #$mub = $vmu[4];
    
    
    printf(OUT "%7.4f %7.4f %7.4f %7.4f %7.4f ",$us,$gs,$rs,$is,$zs);
    printf(OUT "%7.4f %7.4f %7.4f %7.4f %7.4f ",$kp,$V,$R,$I,$wVR);
    printf(OUT "%7.4f %7.4f %7.4f %7.4f %7.4f ",$Z087,$J,$H,$K,$W149);
    printf(OUT "%7.3f %7.3f %7.2f %7.2f %7.2f %7.2f ", $mul,$mub,$Vr,$UU,$VV,$WW);
    printf(OUT "%6.3f %2d %4.2f %5d %5.2f %3d %6.3f %6.3f %8.3f %6.2f ",$Mv,$CL,$Typ,$Teff,$logg,$Age,$Mass,$Mbol,$Radius,$FEH);
    printf(OUT "%9.5f %9.5f %9.5f %9.5f ",$l,$b,$ra,$dec);
    printf(OUT "%6.3f %7.3f %7.3f %7.3f %6.3f %6.3f\n",$D,$X,$Y,$Z,$A_K,$alpha);
    
}

close(IN);
close(OUT);
