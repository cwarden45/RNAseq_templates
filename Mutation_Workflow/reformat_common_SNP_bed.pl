use warnings;
use strict;

my $inputfile = "mm10_common_dbSNP142_2nt.bed";
my $outputfile = "mm10_common_dbSNP142.bed";

#in file that I tested, 2nd position was the SNP position (regardless of strand)

open(OUT, "> $outputfile") || die("Could not open $outputfile!");
open(INPUTFILE, $inputfile) || die("Could not open $inputfile!");
while (<INPUTFILE>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/\r//g;
		$line =~ s/\n//g;
		my @line_info = split("\t",$line);
		my $chr = $line_info[0];
		my $start = $line_info[1];
		my $end = $line_info[2];
		my $snpID = $line_info[3];
		my $score = $line_info[4];
		my $strand = $line_info[5];
		
		$line_info[1] = $end;
		print OUT join("\t",@line_info),"\n";
	}#end while (<INPUTFILE>)
close(INPUTFILE);
close(OUT);

exit;