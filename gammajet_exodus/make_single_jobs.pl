#!/usr/local/bin/perl

$homedir="/direct/phenix+u/workarea/mjuszkie/gammajet_exodus/";

#@files = `cat dummy.list`;

chomp(@files);

$counter = 0;
#foreach $file (@files)
while ($counter < 5000)
{	
    #@run = split(/dst/,$file);

    #print $run[1]."\n";
    #print $run[1]."\n";

 #   print $run[5]. "\n";
    
 #   @seg = split(/\./,$run[6]);
	    
#    $outfile="exec/run_inc_run7_" . $run[5] . "_" . $6;
    $outfile="exec/run_exodus_single_$counter.csh";

#    chomp($run[6]);

    #chomp($run[1]);
    chomp($file);
	    
    open(OUT,"> $outfile") or die "$outfile\n";
    

    
            print OUT "#!/bin/tcsh\n";
    
	    print OUT "echo ' '\n";
	    print OUT "echo 'START: '`date`\n";
	   # print OUT "source ~/change_setup.csh \n";
	   # print OUT "source ~/setldpath.csh . \n";
	   # print OUT "echo \$LD_LIBRARY_PATH\n";
    
    print OUT "cd $homedir\n"; 
    #print OUT "set filep = \`copyfile.pl  $file mjuszkie\`\n";
    #$filep = `root -b -q -x findfile.C\\\(\\\"$file\\\"\\\) | tail -2 | head -1`;
    $filep = $file;
   # print OUT "echo 'writing file'\n";
   # print OUT "echo $filep \n";

    #$houtfile="/direct/phenix+hp/data10/mjuszkie/exodus/cent0/newII/cent0_etas_".$counter.".root";
    $houtfile="/phenix/scratch/mjuszkie/gammajet.exodus/sharkfins/single_".$counter.".root";
    #$noutfile="/phenix/hp/data10/mjuszkie/run7/pisa_eff/analysis/histofiles/ntuple_". $run[1];    
    #$oscarfile="/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/direct/cent0_all_".$counter.".txt";


    #print OUT "root.exe -b <<EOF\n"; 	

print OUT "hostname \n";
print OUT "/phenix/u/workarea/mjuszkie/gammajet_exodus/exodus<<EOF \n";
print OUT "3 \n";
print OUT "1000000\n";
print OUT "$houtfile\n";
print OUT "88 \n";
print OUT "890.7\n";
print OUT "1.0 \n";
print OUT "4 \n";
print OUT "EOF \n";
print OUT "echo ' ' \n";
print OUT "echo 'END: '`date` \n";
print OUT "echo ' ' \n";

    
    
    $counter = $counter +1;
    
}

`chmod -R a+x exec/`; 
