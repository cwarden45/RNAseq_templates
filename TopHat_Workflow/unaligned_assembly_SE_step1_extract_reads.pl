use warnings;
use strict;

my %test_samples = ("",1,"",1);

my $alignment_folder = "/path/to/TopHat_Alignment";
my $reads_folder = "Downsampled_Reads";

opendir DH, $alignment_folder or die "Failed to open $alignment_folder: $!";
my @files= readdir(DH);

foreach my $file (@files){
	my $unaligned_bam = "$alignment_folder/$file/unmapped.bam";
	if((-f $unaligned_bam)&&(exists($test_samples{$file}))){
		print "$file\n";
		
		my $unaligned_FQ = "$reads_folder/$file\_R1.fastq";
		
		my $temp_sam = "$file.bam";
		my $command = "samtools view -h $unaligned_bam | head -n 1000000 > $temp_sam";
		system($command);

		$command = "samtools fastq $temp_sam > $unaligned_FQ";
		system($command);
		
		$command = "rm $temp_sam";
		system($command);
	}
}#end foreach my $file (@files)

exit;