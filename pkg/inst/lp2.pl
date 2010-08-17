#!/usr/bin/perl -w

#Ronnies code to transform Grid output to outbred.qtl input

$filein = "Gout.coe";

open(IN3,$filein); #change to get arguments

@linegrab = ();
@whichline = ();
$linenum = 0;
$fam = 0;
$chromnum = 0;

while (<IN3>){
	chomp;
	@linegrab = split;
	if ($linegrab[0] eq "Family"){
		$whichline[$linenum] = $fam;
		$linenum = $linenum+1;
	}
	if ($linenum > 0){#find the end only when there is something allreay found
		if ($linegrab[0] eq "Chr"){
			$whichline[$linenum] = $fam-1;
			$linenum++
		}
	}
	if ($linegrab[0] eq "#F2"){#end of the file
		$whichline[$linenum] = $fam-1;
		$linenum++
	}
	$fam++;
}

$chromnum = (scalar(@whichline))/2;
print "Data found for ",$chromnum," chromosomes\n";



#tester array finding
#$max = $chromnum;
#print  $max, " array length \n";
#for ($i=0; $i < ($max); $i++) {
#	print   "Item number ",$i, " at line ", $whichline[$i], "\n";
#}

close IN3;	
open(IN3,$filein);

$counter=0;
$counterB = 0;
$counterC = 1;
$setter=0;
@arr =();

#for ($countx = 0; $countx < (1+scalar(@whichline)); $countx++) {
#	print $countx, "\t",$whichline[$countx] ,"\n";
#}


while (<IN3>){
	chomp;
	@arr = split;
	if ($setter==1){
		if (($arr[3]+$arr[4]+$arr[5]+$arr[6]) < 0.98){
			#print "sum ", ($arr[3]+$arr[4]+$arr[5]+$arr[6]),"\n";
			print OUT3 $arr[1], "\t",0.25, "\t",0.5, "\t",0.25,"\n";
		}else{
			print OUT3 $arr[1], "\t",$arr[3], "\t",$arr[4]+$arr[5], "\t",$arr[6],"\n";
		}
	}
	
	#print $whichline[$counterB],"\t"; 

	if ($counterB < scalar(@whichline))	{
		
	if ($counter == $whichline[$counterB]){
		#print  "Out here\n";
		#print $setter , "\n";
		if ($setter==1){
			#	print  "To  FALSE\n";
			$setter = 0;
			#close and open files
			close OUT3;
		} else {
			#print  "To  TRUE\n";
			$setter = 1;
			print  "Converting chromosome ",$counterC,"\n";
			open(OUT3,">p_output_chrom_".$counterC.".txt"); 
			$counterC++;
			
		}
		$counterB++;		
			}
		
	}


	
	$counter++;
}

close IN3;
