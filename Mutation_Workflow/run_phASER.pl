use warnings;
use strict;

my %finished_samples;

my $alignment_folder = "path/to/Alignment_Folder";
my $ref = "/path/to/ref.fa";

my $vcf = "joint_variant_calls.vcf";
#my $gene_bed = "TxDb_[genome]_gene.bed";

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
			my $output = "$sampleID\_phASER";
			$command = "python /opt/phaser/phaser/phaser.py --vcf $vcf.gz --sample $sampleID --bam $bam --o $output --baseq 20 --mapq 10 --paired_end 0 --pass_only 0";
			system($command);
			
			#my $haplotype_counts = "$sampleID\_phASER.haplotypic_counts.txt";
			#$output = "$sampleID\_phASER_gene_ae";
			#$command = "python /opt/phaser/phaser_gene_ae/phaser_gene_ae.py --haplotypic_counts $haplotype_counts --features $gene_bed --o $output";
			#system($command);
		}#end if(!exists($finished_samples{$sampleID})
	}#endif(-f ("$alignment_folder/$file") && ($file =~ /.bam$/))
}#end foreach my $file (@files)

exit;
