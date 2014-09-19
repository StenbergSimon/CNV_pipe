#!/usr/bin/env perl 


#
# cnv_pipeline v 0.3
# Simon Stenberg
# Written to find CNV in yeast, would work just as fine in any other organism
#
# This version includes construction of plots for each scaffold
#
# See Help by typing cnv_pipe.pl -h or --help
#
# this version has a probably better way of calculating normalization.


use warnings;
use strict;
use Cwd;

my @names;
my @mediancov;
my @lengths;
my $line;
my @windowmedians;
my $optionsref =  options(\@ARGV);
my @options = @{$optionsref};
my $bam = $options[0];
my $windowsize = $options[1];
my $cutoff = $options[2];
my $regexp = $options[3];
my $report = $options[4];
my $mapq = $options[5];
my $reference = $options[6];
my $path = getcwd;
my $normalizer = 0;
my $splitbamref;
my $increment = $options[7];
my $samplecovsum = 0;
my $referencecovsum = 0;

our $usedbam = trimfilename($bam);


my $splitbam = $bam;
$splitbam =~ s/.bam$//;
unless($reference eq 0){
	$splitbamref = $reference;
	$splitbamref =~ s/.bam$//;
	our $usedref = trimfilename($reference);
}
# Get names from samtools store in @names

open (ST, "samtools view -H $bam|") || die ($!);
while( $line = <ST>){
 
	if ($line =~ m/$regexp/){
		push (@names, $2);
	}
}
close (ST);

# Make temp dirs and split bam file

if($report ==1){
	system("mkdir -p reports");
	}
system("mkdir -p cnv_pipe_temps");
system("bamtools split -in $bam -reference");
system("mv $splitbam.REF_* cnv_pipe_temps");
unless($reference eq 0){
	system("bamtools split -in $reference -reference");
	system("mkdir -p cnv_pipe_temps/reference");
	system("mv $splitbamref.REF_* cnv_pipe_temps/reference");
	}

# Get lengths store in @lengths

open (ST, "samtools view -H $bam|") || die ($!);
while( $line = <ST>){

	if ($line =~ m/(LN:)(\d+)/){
		push (@lengths, $2);
	}
}

close(ST);

if ($reference eq 0){
	#Median Coverage for each name

	foreach my $name(@names){
		my @templengths;
		open (ST, "samtools view -bq $mapq $bam $name | genomeCoverageBed -d -ibam stdin|") || die($!);
			while ($line = <ST>){
				if($line =~ m/($name)(\s+)(\d+)(\s+)(\d+)/){
					push (@templengths, $5);
					}
			}
		
		my $median = median(\@templengths);
		#$median = log_2($median);
		push (@mediancov, $median);
		close (ST);

	}


	my $namecounter = 0;

	foreach my $name(@names){
		my @temparray;
		my $quote;
		my $entry = " ";
		my $startpos = 0;
		my $endpos = $windowsize;
		my $tempsplitbamname = trimfilename($splitbam) . '.REF_' . $name . '.bam';
		if ($report == 1){
				 open(REPORT, ">reports/$name");
				}
		open(ST, "samtools view -uq $mapq cnv_pipe_temps/$tempsplitbamname | genomeCoverageBed -d -ibam stdin|") || die($!);
		while ($line = <ST>){
			if($line =~ m/($name)(\s+)(\d+)(\s+)(\d+)/){
									push (@temparray, $5);
									}
				}
		until(!defined($temparray[($windowsize-1)])){
			my $i = 0;
			my @windowarray = 0;
			until ($i==($windowsize-1)){
				my $temp = shift(@temparray);
				push(@windowarray, $temp);
				++$i;
				}
			my $tempavg = average(\@windowarray);
			#$tempavg = log_2($tempavg);
			my $diff = (log_2($tempavg) - log_2($mediancov[$namecounter]));
			if ($mediancov[$namecounter] > 0) {
				$quote =( $tempavg / $mediancov[$namecounter]);
				}
				else { $quote = $tempavg;
					 }
			$quote = log_2($quote);
			if ($quote > $cutoff){
				$entry =  $entry . $name . "," . $quote . "," . log_2($mediancov[$namecounter]) . "," . $diff . "," . $startpos . "," . $endpos . ",";

				}
			if ($report ==1){
					
					print REPORT $quote, "\t", $startpos, "\t", $endpos;
					if($quote > $cutoff){
                                                print REPORT "\t", 1, "\n";
                                                }
                                        elsif($quote < ((-1)*$cutoff)){
                                                print REPORT "\t", -1, "\n";
                                                }
                                        else{
                                                print REPORT "\t", 0, "\n";
                                                }
                                        }	
			
			$startpos = $startpos + $increment;
					$endpos = $endpos + $increment;


			}
		$windowmedians[$namecounter] = $entry;
		++$namecounter;
		if($report ==1){
				close(REPORT);
				}

	}
}

else{
	my @refarray;
	my @samplearray;
	my $refreadcount;
	my $samplereadcount;
	open (ST, "samtools view -bq $mapq $bam | genomeCoverageBed -d -ibam stdin|") || die($!);
			while ($line = <ST>){
				if($line =~ m/(\w+)(\s+)(\d+)(\s+)(\d+)/){
					my $push = $5 + 1;
					push (@samplearray, $5);
					}
			}
	open (ST, "samtools view -bq $mapq $reference | genomeCoverageBed -d -ibam stdin|") || die($!);
		while ($line = <ST>){
				if($line =~ m/(\w+)(\s+)(\d+)(\s+)(\d+)/){
					my $push = $5 + 1;
					push (@refarray, $5);
					}
			}
	#Normaliation 1:
	#	my $refmed = median(\@refarray);
	#	my $samplemed = median(\@samplearray);
	#	if ($refmed == 0){
	#		$normalizer = 0;
	#		print "\nWARNING: REFERENCE MEDIAN = 0!\n";
	#	}
	#	else{
	#		$normalizer = $samplemed / $refmed;
	#	}

	#Normalization 2:
	#open (ST, "samtools view -q $mapq -c $reference|") || die($!);
	#	while ($line = <ST>){
	#				$refreadcount = $line;
	#				$refreadcount =~ s/\n//;
	#				$refreadcount =~s/\s//;
	#	}
	#
	#open (ST, "samtools view -q $mapq -c $bam|") || die($!);
        #        while ($line = <ST>){
        #                                $samplereadcount = $line;
        #                                $samplereadcount =~ s/\n//;
        #                                $samplereadcount =~s/\s//;
        #        }
	#$normalizer = $refreadcount / $samplereadcount;
	
	#Normalization 3:
	
	foreach (@samplearray){
		unless($_ == 0){
			$samplecovsum = $samplecovsum + $_;
			}
		} 
	foreach (@refarray){
		unless($_ == 0){
			$referencecovsum = $referencecovsum + $_;
			}
		}
	$normalizer = $referencecovsum / $samplecovsum;		


	my $namecounter = 0;

	foreach my $name(@names){
		my @temparraysample;
		my @temparrayreference;
		my $quote;
		my $entry = " ";
		my $startpos = 0;
		my $endpos = $windowsize;
		my $tempsplitbamname = trimfilename($splitbam) . '.REF_' . $name . '.bam';
		my $tempsplitbamnameref = trimfilename($splitbamref) . '.REF_' . $name . '.bam';
		if ($report == 1){
				 open(REPORT, ">reports/$name");
				}
		open(ST, "samtools view -uq $mapq cnv_pipe_temps/$tempsplitbamname | genomeCoverageBed -d -ibam stdin|") || die($!);
		while ($line = <ST>){
			if($line =~ m/($name)(\s+)(\d+)(\s+)(\d+)/){
									push (@temparraysample, $5);
									}
				}
		
		open(ST,"samtools view -uq $mapq cnv_pipe_temps/reference/$tempsplitbamnameref | genomeCoverageBed -d -ibam stdin|") ||die($1);
		while($line = <ST>){
			if($line =~ m/($name)(\s+)(\d+)(\s+)(\d+)/){
									push (@temparrayreference, $5);
									}
				}
				
		until(!defined($temparraysample[($windowsize-1)])){
			my $i = 0;
			my $i2 = 0;
			my @windowarraysample = 0;
			my @windowarrayreference = 0;
			until ($i==($windowsize-1)){
				my $temp = shift(@temparraysample);
				$temp = $temp * $normalizer;
				push(@windowarraysample, $temp);
				++$i;
				}
			until ($i2==($windowsize-1)){
				my $temp = shift(@temparrayreference);
				push(@windowarrayreference, $temp);
				++$i2;
				}
			my $tempavgsample = average(\@windowarraysample);
			my $tempavgreference = average(\@windowarrayreference);
			if (($tempavgsample > 0) && ($tempavgreference > 0)) {
				$quote =( $tempavgsample / $tempavgreference);
				}
				elsif(($tempavgsample == 0) && ($tempavgreference > 0 )){
					$quote = 0.1;
					}
				else { $quote = $tempavgsample;
					 }
			$quote = log_2($quote);
			if (($quote > $cutoff) || $quote < (((-1)*$cutoff))){
				$entry =  $entry . $name . "," . $quote . "," . log_2($tempavgreference) . "," . log_2($tempavgsample) . "," . $startpos . "," . $endpos . ",";

				}
			if ($report ==1){
					
					print REPORT $quote, "\t", $startpos, "\t", $endpos;
					if($quote > $cutoff){
						print REPORT "\t", 1, "\n";
						}
					elsif($quote < ((-1)*$cutoff)){
						print REPORT "\t", -1, "\n";
						}
					else{
						print REPORT "\t", 0, "\n";
						}
					}	
			
			$startpos = $startpos + $increment;
					$endpos = $endpos + $increment;


			}
		$windowmedians[$namecounter] = $entry;
		++$namecounter;
		if($report ==1){
				close(REPORT);
				}

		}
	}


#Output:

if ($reference eq 0){
	output_results(\@windowmedians);	
	}
else{
	output_results_with_ref(\@windowmedians);
	}
	
#Plotting

if ($report ==1){
		foreach(@names){ 
				system("Rscript ~/bin/plot_cnv.r reports/$_");
			}
		system("rm $path/cnv_pipe_temps/output");
		system("mkdir -p reports/plots");
                system("mv reports/*.pdf reports/plots");

	}

# Remove temp folder with files

system("rm -r cnv_pipe_temps");


##########################################################################################
#####Subroutines##########################################################################
##########################################################################################


sub help {
	die( "\n -------------------- cnv_pipeline v:0.3 ------------------------------
	\nDependencies: Bedtools, samtools, bamtools.
	\nUsage: cnv_pipe.pl sortedbamfile.bam [Options] > output.csv
	\nBAM-file(s) _must_ be sorted and indexed.
	\nIncreasing the windowsize, increases the runtime by equal proportions.
	\nIt's highly recommended to consider GC-coverage-bias in your sample and reference
	\n
	\nOptions:
	\n-w [int]\t # Sliding windowsize. Mind the total depth and significance level. Also determines sensetivity. Default: 500
	\n-co [int]\t # Determines how much greater in log2 scale the regions has to be in depth to be called. Default: 0.4
	\n-rall\t\t #Prints out all log2 ratios regardless if cutoff is met. Results in a tab-separated file for each chr/scf
	\n-incr [int]\t # Window increment size, default is = windowsize (non-overlapping windows)
	\n-mapq [int]\t #  MAPQ threshold for calculating coverage, Defualt = 0 (will include reads mapping to multiple places)
	\n-r reference.bam # Run the program using comparison to reference alignment (control) vs sample alignment
	\n-regex [Reg-ex]\t # !!!NOT FINISHED, May not work!!!Perl-style regular expression (excluding //) to match the names in your bam file.
	\n\t\t #Regular expression must match the name in group 2, so group 1 should always be (SN:). ", 'Default: (SN:)(\D+\d+)', "
	\n-----------------------------------------------------------------------
	\n");
	}
	
sub median{

	my @array_unsorted = @{$_[0]};

	my @array = sort {$a <=> $b} @array_unsorted;
	my $pos = ((scalar(@array)/2) - 1);
	
	unless(defined($array[$pos]) && defined($array[($pos+1)])){
	
		return "0";
		}
	
	else{
	
		 if( (scalar(@array)) % 2){

			 
			 my $median = $array[($pos+0.5)];
			 return $median;
	
	
		 } else {
	
			 my $pos = ((scalar(@array)/2) - 1);
	
			 my $tot = ($array[$pos] + $array[($pos+1)]);
			 my $median = $tot/2;
	
			 return $median;
	
			 }
	}
}

sub log_2{

  my $num = shift;
  my $base = 2;
  if ($num > 0){
  	return log($num)/log($base);
	}
	else{
		return 0;
		}
}

sub average{

	my @array = @{$_[0]};
	my $sum = 0;
	my $mean = 0;
	
	foreach (@array){
	
		$sum = $sum + $_ ; 
	
	}
	
	unless(scalar(@array) == 0){
	
		$mean = $sum / (scalar(@array));
		return $mean;
		}	
	else{
		

	return 0;
	}
}


sub output_results{
	# Use: output_results(@windowmedians)
	my @entries = @{$_[0]};
	my $i = 0;
	
	print "#CNV_Analysis\n#Sample:$main::usedbam\n#Windowsize:$windowsize\n#Cutoff:$cutoff\n#Increment:$increment\n#Mode=One-sample\n";
	print "#Chrom/Contig\tLog2Ratio\tContigMedian\tDifference\tWindowStart\tWindowStop\n";
	
	foreach (@entries){
		if($_ ne " "){
			$_ =~ s/\s//;
			
			my @array = split("," , $_);
			my $i2 = 0;	
			while(($i2 + 6) <= scalar(@array)){
					
					print $array[$i2];
					print "\t";
					print $array[($i2+1)];
					print "\t";
					print $array[($i2+2)];
					print "\t";
					print $array[($i2+3)];
					print "\t";
					print $array[($i2+4)];
					print "\t";
					print $array[($i2+5)];
					print "\n";
					
					$i2 = $i2 + 6;
				}
			}
		}
	++$i;
}
sub output_results_with_ref{
	# Use: output_results_with_ref(@windowmedians)
	my @entries = @{$_[0]};
	my $i = 0;
	
	print "#CNV_Analysis\n#Sample:$main::usedbam\n#Reference:$main::usedref\n#Windowsize:$windowsize\n#Cutoff:$cutoff\n#Increment:$increment\n#Mode=Control\n";
	print "#Chrom/Contig\tLog2_Ratio\tReference_Average\tSample_Average\tWindowStart\tWindowStop\n";
	
	foreach (@entries){
		if($_ ne " "){
			$_ =~ s/\s//;
			
			my @array = split("," , $_);
			my $i2 = 0;	
			while(($i2 + 5) <= scalar(@array)){
					
					print $array[$i2];
					print "\t";
					print $array[($i2+1)];
					print "\t";
					print $array[($i2+2)];
					print "\t";
					print $array[($i2+3)];
					print "\t";
					print $array[($i2+4)];
					print "\t";
					print $array[($i2+5)];
					print "\n";
					
					$i2 = $i2 + 6;
				}
			}
		}
	++$i;
}



sub options{ 
	my $windowsize;
	my $bam;
	my $cutoff;
	my $regexp;
	my $i2 = 0;
	my @options = @{$_[0]};
	my @sender;
	my $report;
	my $mapq;
	my $reference;
	my $increment;
	
	foreach(@options){
	

		if($_ eq "-w"){
						$windowsize = $options[($i2+1)];
								if($windowsize =~ m/\D/){
								die("\nError: Illegal windowsize, has to be numerical\n", help());
								}
				
						}
		elsif($_ eq "-co"){
						 $cutoff = $options[$i2+1];
						unless (defined $cutoff && $cutoff =~ m/\d+/ ){
							die("Error: Need cutoff value if option is set\n", help());
							}
						  }
					
		
		elsif(($_ eq "-h") or ($_ eq "--help")){
						
						die(help());
					
						}
		elsif($_ eq "-regex"){
						$regexp = $options[$i2+1];
						#Group 2 has to match name, first group should always be (SN:).
						}
		elsif($_ eq "-rall"){
						$report = 1;
					}
		elsif($_ eq "-mapq"){
					$mapq = $options[$i2+1];
				
				}
		elsif($_ eq "-r"){
				$reference  = $options[$i2+1];	
				}
		elsif($_ eq "-incr"){
				$increment = $options[$i2+1];
				}
			
		$i2++;		
		
	
		}

	# First argument should be file name - load into $file

	$bam = $ARGV[0];

	# No file? Die with helpmsg

	unless (defined $bam){
		die(help());
		}
	
	unless (defined $windowsize){
		$windowsize = 500;	
		}
	unless (defined $cutoff){
		$cutoff = 0.4;
		}
	unless (defined $regexp){
		$regexp = '(SN:)([^\s]+)';
		}
	unless (defined $report){
			$report = 0;
		}
	unless(defined $mapq){
		$mapq = 0;
		}
	unless(defined $increment){
		$increment = $windowsize;
		}
	

	# Wrong fileformat? die.

	unless( $bam =~ m/\.bam$/){
		die ( "\nNeed .bam as input", help() );
		}
	if (defined $reference){
		unless ($reference =~ m/\.bam$/){
			die ("\nReference has to be in .bam format")
			}
		}

	unless(defined $reference){
                $reference = 0;
                }



	# Check that all arguments submitted are valid

	foreach(@ARGV){
		unless($_ =~ m/-h/ || $_ =~ m/--help/ || $_ =~ m/-regex/ || $_ =~ m/-co/ || $_ =~ m/-w/ || $_ =~ m/$bam/ || $_ =~ m/\Q$cutoff/ || $_ =~ m/\Q$windowsize/ || $_ =~ m/\Q$regexp/ || $_ =~ m/-rall/ || $_ =~ m/$report/ || $_ =~ m/$mapq/ || $_ =~ m/-mapq/ || $_ =~ m/-r/ || $_ =~ m/$reference/ || $_ =~ m/-incr/ || $_ =~ m/$increment/){
			print "\n Error: Invalid arguments\n ";
			die (help());
			}
		}
	
	push(@sender, $bam);	
	push(@sender, $windowsize);
	push(@sender, $cutoff);
	push(@sender, $regexp);
	push(@sender, $report);
	push(@sender, $mapq);
	push(@sender, $reference);
	push(@sender, $increment);

	return (\@sender);

}	

	
sub trimfilename{
	
	my $file = $_[0];

	my @array = split("/",$file);

	return ($array[(scalar(@array)-1)]);
}




