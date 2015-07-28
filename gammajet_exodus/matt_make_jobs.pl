#!/usr/local/bin/perl

$homedir="/direct/phenix+u/workarea/mjuszkie/gammajet_exodus/";

#@files = `cat dummy.list`;

#chomp(@files);

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

#HERE:The exectuable names: 
    $outfile="exec/run_exodus_$counter.csh";

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

    $filep = $file;


##HERE: Change houtfile where you want the root files to go; oscarfile doesn't actually get written out at all so path shouldn't matter
    $houtfile="/direct/phenix+hp/data10/mjuszkie/exodus/cent0/new/cent0_etas_".$counter.".root";
    $oscarfile="/phenix/scratch/mjuszkie/gammajet.exodus/Takao/Cent0/direct/cent0_all_".$counter.".txt"; 


    #print OUT "root.exe -b <<EOF\n"; 	

print OUT "hostname \n";
print OUT "/phenix/u/workarea/mjuszkie/gammajet_exodus/exodus<<EOF \n";
print OUT "6 \n";
print OUT "50000 \n";
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
