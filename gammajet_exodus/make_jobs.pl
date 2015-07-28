#!/usr/local/bin/perl

$homedir="/phenix/hhj/ahanks/gammajet_exodus";

$counter = 0;
while ($counter < 4000)
{	
    $outfile="exec/run_exodusII_$counter.csh";

    open(OUT,"> $outfile") or die "$outfile\n";
        
	print OUT "#!/bin/tcsh\n";
    
	print OUT "echo ' '\n";
    	print OUT "echo 'START: '`date`\n";
    
    	print OUT "cd $homedir\n"; 

    $houtfile="$homedir/fe_sharks/cent0_II_".$counter.".root";
    $oscarfile="$homedir/fe_sharks/cent0_all_".$counter.".txt";

	print OUT "hostname \n";
	print OUT "$homedir/exodus<<EOF \n";
	print OUT "6 \n";
	print OUT "25000 \n";
	print OUT "$houtfile\n";
	print OUT "$oscarfile\n";
	print OUT "88 \n";
	print OUT "6 \n";
	print OUT "EOF \n";
	print OUT "echo ' ' \n";
	print OUT "echo 'END: '`date` \n";
	print OUT "echo ' ' \n";
    
    $counter = $counter +1;
    
}

`chmod -R a+x exec/`; 
