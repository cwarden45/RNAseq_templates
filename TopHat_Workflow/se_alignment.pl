#Code written by Charles Warden (cwarden@coh.org, x60233)

use warnings;
use strict;

my $root = "/path/to/folder";
my $ref = "/path/to/ref";
my $threads = 4;

my $SE_reads_folder = "$root/Reads";

my $alignment_folder = "$root/TopHat_Alignment";

my @finished_files = ();
RNA_Seq_workflow_se_no_gz(\@finished_files, $SE_reads_folder, $alignment_folder, $threads, $ref);
exit;
	
sub RNA_Seq_workflow_se_no_gz
	{
		my ($arr_ref, $reads, $outputfolder, $threads, $ref)=@_;
		my @include_samples = @$arr_ref;
		
		#standard stuff
		
		opendir DH, $reads or die "Failed to open $reads: $!";
		my @files = readdir(DH);
		foreach my $file (@files)
			{
				#print "$file\n";
				if(-f ("$reads/$file") && ($file =~ /_L\d{3}_R1_001.fastq/))
					{
						my ($sampleID)= ($file =~ /(.*)_\w{6}_L\d{3}_R1_001.fastq/);
						print "Working on $sampleID\n";
										
						my $sample_folder = "$sampleID";
									
						my $new_flag = 1;
							
						foreach my $old_sample (@include_samples)
							{
								if ($old_sample eq $sample_folder)
										{
											$new_flag = 0;
										}
							}#end foreach my $old_sample (@exclude_samples)
								
						if($new_flag)
							{
								my $output_subfolder="$outputfolder/$sample_folder";
								mkdir($output_subfolder);
									
								my $read1 = "$reads/$file";
								
								print "\n\nAlign via TopHat\n\n";
								my $command = "tophat -o $output_subfolder -p $threads --no-coverage-search $ref $read1";
								system($command);
								
								my $bam_file = "$output_subfolder/accepted_hits.bam";																			
								my $rename_bam = "$outputfolder/$sampleID.bam";
								
								print "\n\nCreate Sorted BAM File\n\n";
								$command = "/opt/samtools-1.3/bin/samtools sort $bam_file -o $rename_bam";
								system($command);

								$command = "rm $bam_file";
								system($command);

								print "\n\nIndexing BAM File\n\n";
								$command = "/opt/samtools-1.3/bin/samtools index $rename_bam";
								system($command);
							}#end if if($new_flag)
					}#end if(-f ("$inputfolder/$file") && ($file =~ /_R1_/))
			}#end foreach my $file (@files)
		closedir(DH);
	}#end DNA_target_workflow