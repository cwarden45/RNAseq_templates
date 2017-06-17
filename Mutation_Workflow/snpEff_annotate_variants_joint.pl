use warnings;
use strict;

my $jointVCF = "joint_variant_calls.GATK.HC.best.practices.filtered.sansRNAedit.vcf";
my $resultFolder = "Annotated_Variants/GATK_Joint/";

my $parameter_file = "parameters.txt";
my $java_mem = "";
my $snpEff_path = "";
my $genome = "";
my $dbSNP_bed = "";
my $sample_file = "";

open(INPUTFILE, $parameter_file) || die("Could not open $parameter_file!");
while (<INPUTFILE>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		my @line_info = split("\t",$line);
		my $param = $line_info[0];
		my $value = $line_info[1];
		
		if ($param eq "Java_Mem"){
			$java_mem = $value;
		}
			
		if ($param eq "SnpEff_Path"){
			$snpEff_path = $value;
		}

		if ($param eq "genome"){
			$genome = $value;
		}	
		
		if ($param eq "dbSNP_BED"){
			$dbSNP_bed = $value;
		}
		
		if ($param eq "sample_description_file"){
			$sample_file = $value;
		}	
		
	}#end while (<INPUTFILE>)
close(INPUTFILE);

if (($java_mem eq "")||($java_mem eq "[Required]")){
	die("Need to enter a value for 'Java_Mem'!")
}

if (($snpEff_path eq "")||($snpEff_path eq "[Required]")){
	die("Need to enter a value for 'SnpEff_Path'!")
}

if (($genome eq "")||($genome eq "[Required]")){
	die("Need to enter a value for 'genome'!")
}

if (($dbSNP_bed eq "")||($dbSNP_bed eq "[Required]")){
	die("Need to enter a value for 'dbSNP_BED'!")
}

if (($sample_file eq "")||($sample_file eq "[Required]")){
	die("Need to enter a value for 'sample_description_file'!")
}

my $ORegAnno_bed = "$snpEff_path/$genome\_ORegAnno.bed";

my $vcf_out = "$resultFolder/snpEff_annotations.vcf";
my $html_out = "$resultFolder/snpEff_annotations.html";
my $csv_out = "$resultFolder/snpEff_statistics.csv";

my $command = "java -jar -Xmx$java_mem $snpEff_path/snpEff/snpEff.jar $genome $jointVCF -interval $ORegAnno_bed -interval $dbSNP_bed -csvStats $csv_out -stats $html_out > $vcf_out";
#system($command);

#create short ID mapping
open(META, $sample_file) || die("Could not open $sample_file!");

my $line_count = 0;
my $long_index = -1;
my $short_index = -1;
my %sample_hash;

while (<META>){
	my $line = $_;
	chomp $line;
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	$line_count++;
	
	my @line_info = split("\t",$line);
	
	if ($line_count == 1){
		for (my $i=0; $i < scalar(@line_info); $i++){
			if($line_info[$i] eq "sampleID"){
				$long_index = $i;
			}elsif($line_info[$i] eq "userID"){
				$short_index = $i;
			}
		}#end for (my $i=0; $i < scalar(@line_info); $i++)
	}elsif(($long_index != -1)&&($short_index != -1)){
		$sample_hash{$line_info[$long_index]}=$line_info[$short_index]
	}else{
		print "Your sample description file must have 'sampleID' (based upon file names) and 'userID' (name desired in table, which should not start with number\n";
		print "Also, you probably want to comment out line running snpEff\n";
		exit;
	}#end else
		
}#end while (<META>)

close(META);

#reformat .vcf to tab-delimited text with values of interest
my $revised_table = "$resultFolder/snpEff_annotated_genotypes.txt";
open(OUT, "> $revised_table") || die("Could not open $revised_table!");

open(IN, $vcf_out) || die("Could not open $vcf_out!");
while (<IN>){
	my $line = $_;
	chomp $line;
	$line =~ s/\r//g;
	$line =~ s/\n//g;
	if(!($line =~ /^##/)){
		my @line_info = split("\t",$line);
			
		if($line =~ /^#/){
			my $header = "Chr\tPos\tRef\tAlt";
			for (my $i=9; $i < scalar(@line_info); $i++){
				if(exists($sample_hash{$line_info[$i]})){
					$header = "$header\t".$sample_hash{$line_info[$i]}.".genotype";
				}else{
					print "Could not find ".$line_info[$i]." in sample mapping.  Please revise sample description file\n";
					print "Also, you probably want to comment out line running snpEff\n";
					exit;
				}#end else
			}#end for (my $i=9; $i < scalar(@line_info); $i++)
			$header = "$header\tAnnotated.Nuc\tLocation\tRefGene\tOverall.Flag\tPred.Effect\tGene.AA\tGene.Nuc\tAnn.Dist\tORegAnno\tCommon.dbSNP\n";
			print OUT $header;
		}else{
			my $chr = $line_info[0];
			my $pos = $line_info[1];
			my $ref = $line_info[3];
			my $var = $line_info[4];
			my $text = "$chr\t$pos\t$ref\t$var";
			for (my $i=9; $i < scalar(@line_info); $i++){
				my @geno_info = split(":",$line_info[$i]);
				$text = "$text\t".$geno_info[0];
			}#end for (my $i=9; $i < scalar(@line_info); $i++)
			
			my @var_info = split(";",$line_info[7]);
			foreach my $info (@var_info){
				if($info =~ /^ANN=/){
					$info =~ s/^ANN=//;
					my @separate_ann = split(",",$info);
					my @ann_info = split("\\|",$info);
					my $nuc = "";
					my $location = "";
					my $impact = "";
					my $gene = "";
					my $dist = "";
					my $gene_nuc = "";
					my $gene_AA = "";
					my $ORegAnno = "NA";
					my $dbSNP = "NA";
					
					for (my $i=0; $i < scalar(@separate_ann); $i++){
						my $snpEff_ann = $separate_ann[$i];
						my @ann_info = split("\\|",$snpEff_ann);
						my $temp_nuc = $ann_info[0];
						my $temp_location = $ann_info[1];
						my $temp_impact = "UNKNOWN";
						my $temp_gene = "NA";
						my $temp_dist = "NA";
						my $temp_gene_nuc = "NA";
						my $temp_gene_AA = "NA";
						
						if ($snpEff_ann =~ /ORegAnno/){
							if ($ORegAnno eq "NA"){
								$ORegAnno=$ann_info[6];
							}else{
								#print "Multiple ORegAnno Annotations for variant...\n";
								$ORegAnno=$ORegAnno.",".$ann_info[6];
							}
						}elsif($snpEff_ann =~ /dbSNP/){
							if ($dbSNP eq "NA"){
								$dbSNP=$ann_info[6];
							}else{
								#print "Multiple Common dbSNP Annotations for variant...\n";
								$dbSNP=$dbSNP.",".$ann_info[6];
							}
						}else{
							if (scalar(@ann_info) > 3){
								$temp_impact = $ann_info[2];
								$temp_gene = $ann_info[3];
								$temp_gene_nuc = $ann_info[6].":".$ann_info[9];
								if (scalar(@ann_info) > 10){
									$temp_gene_AA = $ann_info[10];
									if(scalar(@ann_info) >=15){
										$temp_dist = $ann_info[14];
									}#end if(scalar(@ann_info) >=14)
								}#end if (scalar(@ann_info) > 9)
							}#end if (scalar(@ann_info) > 3)
							
							if ($i == 0){
								$nuc=$temp_nuc;
								$location=$temp_location;
								$impact=$temp_impact;
								$gene = $temp_gene;
								$dist = $temp_dist;
								$gene_nuc=$temp_gene_nuc;
								$gene_AA = $temp_gene_AA;
							}else{
								$nuc="$nuc,$temp_nuc";
								$location="$location,$temp_location";
								$impact="$impact,$temp_impact";
								$gene = "$gene,$temp_gene";
								$dist = "$dist,$temp_dist";
								$gene_nuc="$gene_nuc,$temp_gene_nuc";
								$gene_AA ="$gene_AA,$temp_gene_AA";
							}#end else
						}#end else
					}#end for (my $i=0; $i < scalar(@separate_ann); $i++)

					$impact =~ s/MODIFIER/UNKNOWN/g;
					my $overall_flag = "UNKNOWN";
					if($dbSNP ne "NA"){
						$overall_flag="COMMON";
					}elsif($impact =~ /HIGH/){
						$overall_flag="HIGH";
					}elsif($impact =~ /MODERATE/){
						$overall_flag="MODERATE";
					}elsif($impact =~ /LOW/){
						$overall_flag="LOW";
					}elsif($ORegAnno ne "NA"){
						$overall_flag="ORegAnno";
					}else{
						my $filtered_impact = $impact;
						$filtered_impact =~ s/UNKNOWN//g;
						$filtered_impact =~ s/,//g;
						if ($filtered_impact ne ""){
							print "Modify code to give assignment for $filtered_impact / $impact\n";
							exit;
						}
					}#end else
					$text = "$text\t$nuc\t$location\t$gene\t$overall_flag\t$impact\t$gene_AA\t$gene_nuc\t$dist\t$ORegAnno\t$dbSNP\n";
					print OUT $text;
				}#end if($info =~ /^ANN=/)
			}#end foreach my $info (@var_info)
		}
	}#end if(!($line =~ /^##"/))
		
}#end while (<INPUTFILE>)
close(IN);
close(OUT);

exit;
