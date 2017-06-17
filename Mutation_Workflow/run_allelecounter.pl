use warnings;
use strict;

#Download from https://github.com/secastel/allelecounter

my %finished_samples;

my $alignment_folder = "path/to/Alignment_Folder";
my $ref = "/path/to/ref.fa";

my $vcf = "joint_variant_calls.vcf";

my $command = "bgzip $vcf";
system($command);

$command = "tabix -p vcf $vcf.gz";
system($command);

opendir DH, $alignment_folder or die "Failed to open $alignment_folder: $!";
my @files = readdir(DH);
foreach my $file (@files){
	if(-f ("$alignment_folder/$file") && ($file =~ /.bam$/)){
		my ($sampleID)= ($file =~ /(.*).bam/);
		print "Working on $sampleID\n";
		
		if(!exists($finished_samples{$sampleID})){
			my $bam = "$alignment_folder/$file";
			my $output = "$sampleID\_allele_counter";
			$command = "python /opt/allelecounter/allelecounter.py --vcf $vcf.gz --sample $sampleID --bam $bam --o $output --ref $ref --min_cov 4 --min_baseq 20 --min_mapq 10 --max_depth 10000";
			system($command);
		}#end if(!exists($finished_samples{$sampleID})
	}#endif(-f ("$alignment_folder/$file") && ($file =~ /.bam$/))
}#end foreach my $file (@files)

exit;
