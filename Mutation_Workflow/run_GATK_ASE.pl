use warnings;
use strict;

my %finished_samples;

my $alignment_folder = "path/to/Alignment_Folder";
my $ref = "/path/to/ref.fa";

my $vcf = "joint_variant_calls.vcf";

opendir DH, $alignment_folder or die "Failed to open $alignment_folder: $!";
my @files = readdir(DH);
foreach my $file (@files){
	if(-f ("$alignment_folder/$file") && ($file =~ /.bam$/)){
		my ($sampleID)= ($file =~ /(.*).bam/);
		print "Working on $sampleID\n";
		
		if(!exists($finished_samples{$sampleID})){
			my $bam = "$alignment_folder/$file";
			my $output = "$sampleID\_GATK_ASE.tsv";
			my $command = "java -Xmx16g -jar /opt/GATK-3.7/GenomeAnalysisTK.jar -T ASEReadCounter -R $ref -I $bam -o $output -sites $vcf -U ALLOW_N_CIGAR_READS -minDepth 10 --minMappingQuality 10 --minBaseQuality 20";
			system($command);
		}#end if(!exists($finished_samples{$sampleID})
	}#endif(-f ("$alignment_folder/$file") && ($file =~ /.bam$/))
}#end foreach my $file (@files)

exit;
