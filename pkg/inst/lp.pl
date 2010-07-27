#!/usr/bin/perl -w

# takes file of Carls probabilities and writes as file for each chromosome where each row an individual and the a and d indicator variables are given for each marker in turn
# this is for more than one chromosome and makes a series of files, one per chromosome
# file names are given as command line arguments
# is for F2 but should be ok for backcross too although won't have any d values NB backcross acts similarly to sex chromosome

$Carl_in="cnout";#$ARGV[0]; # output file from cnF2freq
#$pheno_in=$ARGV[1]; # phenotype files for getting sex
#$col_sex=$ARGV[2];# column number for sex in the phenotype file
$sex_chrom= 0;#$ARGV[3]; # which number chromosome is the sex chromosome is zero if sex chromosome not included
$het_sex=$ARGV[4]; # which sex has only one copy of sex chromosome (heterogametic sex)
$ad_values= "p_output";#$ARGV[5]; # main name for output files without .txt ending as are making a series of them

#open (IN2,$pheno_in);
#<IN2>;
#while (<IN2>) {
	#chomp;
	#@in2=split;	
	#$sex{$in2[0]}=$in2[$col_sex-1];	
	#}
#close IN2;

open (IN,$Carl_in);

<IN>;
$f=1;
$p="NA";
while (<IN>) {
	chomp;
	@in=split;
	if ($in[1] && !$in[2]) { # new chromosome	
		if ($f==1) { # very start of file
			$f=0;
		} else {# report for previous chromosome
			close OUT;			
			print "Values calculated for ",$p," positions\n";	
			if ($p-1 != $n) {
				print "Error, the number of positions is not what was expected\n";
			}
		}
		$fs=1; # first individual for chromosome
		$ch=$in[0];
		$n=$in[1];
		open (OUT,">".$ad_values."_chrom_".$ch.".txt");
	} elsif (!$in[0]) { # space between successive individuals	
		if ($fs==1) { # checks if is first individual for that chromosome
			$fs=0;
			#	print OUT "\n";
		} else {
			#print OUT "\n";		
			#		print "Values calculated for ",$p," positions\n";	
		}
	} elsif (!$in[1]) {
		$id=$in[0];
		#print OUT $id;
		$p=0;	
	} else {
		$sum=$in[0]+$in[1]+$in[2]+$in[3];
		$rounded=sprintf("%.3f",$sum); # rounds to 3 digits after decimal point 
		$p++;
		if ($rounded !=1.000) {
			print "Error, for chromosome ",$ch," individual ",$id," does not have a total probability of 1 at position ",$p,"\n";
		}
		unless ($ch==$sex_chrom) {
			$a=$in[0]-$in[3];
			$d=$in[1]+$in[2];
		}
		if ($ch==$sex_chrom) { # for sex chromosome one chromosome is inherited intact from heterogametic parent therefore line origin from one parent is fixed and should only have two columns of non zero probabilities 
			if ($p==1) {
				$het="NA";
				$hom="NA";
			}
			$check=0;
			foreach $t (0..3) {
				if ($in[$t]!=0.000000) {
					$check++;
				}
			}
			if ($check>2) { # if two columns are not zero writes error to screen	
				@order=sort {$a <=> $b} @in[0..3];
				if ($order[1]>0.000002) {
					print "Error individual ",$id," has more than two possible genotypes at position ",$p," on the sex chromosome\n";
					print "The smallest values are: ",$order[0],"\t",$order[1],"\n";
				}
				unless (exists $more_sex{$id}) {
					$more_sex{$id}=$sex{$id};	
				}
			} 
			if ($in[0]!=0.000000) {
				if ($hom eq "NA") {
				      $hom=0; # defines which homozygote column is used
			      	} elsif ($hom==3) {
					print "Error the genotype probabilities for individual ",$id," on the sex chromosome appear to be wrong\n";
				}
				if ($in[3]!=0.000000) { # for sex chromosome shouldn't have both of these columns non zero
					print "Error the genotype probabilities for individual ",$id," at position ",$p," on the sex chromosome appear to be wrong\n";
				}
				if ($in[1]!=0.000000) {
					if ($het eq "NA") {
						$het=1;
					} elsif ($het==2) {
						print "Error the genotype probabilities for individual ",$id," on the sex chromosome appear to be wrong\n";
					}
					if ($in[2]!=0.000000) {
						print "Error the genotype probabilities for individual ",$id," at position ",$p," on the sex chromosome appear to be wrong\n";
					} elsif ($sex{$id}==$het_sex)  { # calculations are different for the heterogametic sex
						$a=($in[0]-$in[1])/2;
						$d=0;
					}
				} elsif ($in[2]!=0.000000) {
					if ($het eq "NA") {
						$het=2;
					} 
					if ($het==1) {
						print "Error the genotype probabilities for individual ",$id," on the sex chromosome appear to be wrong\n";
					} elsif ($sex{$id}==$het_sex) {
						$a=($in[0]-$in[2])/2;
						$d=0;
					}
				} elsif ($sex{$id}==$het_sex) {
					$a=$in[0]/2;
					$d=0;
				}
			} elsif ($in[3]!=0.000000) {
				if ($hom eq "NA") {
					$hom=3;
				} elsif ($hom==0) {
					print "Error the genotype probabilities for individual ",$id," on the sex chromosome appear to be wrong\n";
				}
				if ($in[1]!=0.000000) {
					if ($het eq "NA") {
						$het=1;
					} elsif ($het==2) {
						print "Error the genotype probabilities for individual ",$id," on the sex chromosome appear to be wrong\n";
					}
					if ($in[2]!=0.000000) {
						print "Error the genotype probabilities for individual ",$id," at position ",$p," on the sex chromosome appear to be wrong\n";
					} elsif ($sex{$id}==$het_sex) {
						$a=($in[1]-$in[3])/2;
						$d=0;
					}
				} elsif ($in[2]!=0.000000) {
					if ($het eq "NA") {
						$het=2;
					}
					if ($het==1) {
						print "Error the genotype probabilities for individual ",$id," on the sex chromosome apear to be wrong\n";
					} elsif ($sex{$id}==$het_sex) {
						$a=($in[2]-$in[3])/2;
						$d=0;
					}
				} elsif ($sex{$id}==$het_sex) {
					$a=-$in[3]/2;		
					$d=0;
				}
			} elsif ($in[1]!=0.000000) {
				if ($het eq "NA") {
					$het=1;
				} elsif ($het==2) {
					print "Error the genotype probabilities for individual ",$id," on the sex chromosome appear to be wrong\n";
				} 
				if ($in[2]!=0.000000) {
					print "Error the genotype probabilities for individual ",$id," at position ",$p," on the sex chromosome appear to be wrong\n";
				} elsif ($sex{$id}==$het_sex) {
					$a=-$in[1]/2;
					$d=0;
				}
			} elsif ($in[2]!=0.000000) {
				if ($het eq "NA") {
					$het=2;
				} 
				if ($het==1) {
					print "Error the genotype probabilities for indiviudal ",$id," appear to be wrong for the sex chromosome\n";
				} elsif ($sex{$id}==$het_sex) {
					$a=-$in[2]/2;
					$d=0;
				}
			} else {
			       print "Error individual ",$id," has zero probabilities for all genotypes at position ",$p," on the sex chromosome\n";
		       	}	       
		       	if ($sex{$id}!=$het_sex) {
					$a=$in[0]-$in[3];
					$d=$in[1]+$in[2];
		       	}
	       	}		
		$ab=$in[1]+$in[2];
		print OUT $id, "\t",$in[0],"\t",$ab,"\t",$in[3],"\n";
	}
}
#foreach $u (keys %more_sex) {
#	print "Individual ",$u," sex ",$more_sex{$u}," has more than two genotypes for the Z chromosome\n";
#}
#print "the number of individuals with incorrect genotypes on the sex chromosome is ",scalar(keys %more_sex),"\n";
