use warnings;
use strict;

my $inputfile = "RADAR_mm9.txt";
my $mm9_bed = "mm9_RADAR.bed";
my $bed_for_liftover = "RADAR_mm9_2nt.bed";
my $liftOver = "mm9ToMm10.over.chain";
my $mm10_liftOver = "RADAR_mm10_2nt.bed";
my $unmapped = "RADAR_mm9_mm10_unmapped";
my $mm10_liftOver_reformat = "mm10_RADAR.bed";

my $line_count = 0;

open(BED9, "> $mm9_bed") || die("Could not open $mm9_bed!");
open(BED9b, "> $bed_for_liftover") || die("Could not open $bed_for_liftover!");

open(IN, $inputfile) || die("Could not open $inputfile!");

while (<IN>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/\r//g;
		$line =~ s/\n//g;
				
		$line_count++;

		if ($line_count > 1){
			my @line_info = split("\t",$line);
			my $chr = $line_info[0];
			my $pos = int($line_info[1]);
			my $strand = $line_info[3];
			
			
			print BED9 "$chr\t$pos\t$pos\t.\t0\t$strand\n";
			
			my $pos1;
			my $pos2;
			
			if ($strand eq "+"){
				$pos1 = $pos;
				$pos2 = $pos + 1;
			}elsif($strand eq "-"){
				$pos1 = $pos - 1;
				$pos2 = $pos;
			}else{
				print "$strand not either '+' or '-'\n";
				exit;
			}
			
			print BED9b "$chr\t$pos1\t$pos2\tmm9:$chr:$pos:$strand\t0\t$strand\n";
			
		}#end if ($line_count > 1)
	}#end while (<IN>)
close(IN);
close(BED9);
close(BED9b);

my $command = "liftOver $bed_for_liftover $liftOver $mm10_liftOver $unmapped";
system($command);

open(OUT, "> $mm10_liftOver_reformat") || die("Could not open $mm10_liftOver_reformat!");

open(IN, $mm10_liftOver) || die("Could not open $mm10_liftOver!");

while (<IN>)
	{
		my $line = $_;
		chomp $line;
		$line =~ s/\r//g;
		$line =~ s/\n//g;
				
		my @line_info = split("\t",$line);
		my $chr= $line_info[0];
		my $start = $line_info[1];
		my $stop = $line_info[2];
		my $strand = $line_info[5];
		my $mm9_info = $line_info[3];
		
		my ($mm9_strand) = ($mm9_info =~ /:(\S$)/);
		
		if ($strand eq $mm9_strand){
			if ($mm9_strand eq "+"){
				my $pos = $start;
				print OUT "$chr\t$pos\t$pos\t.\t0\t$strand\n";
			}else{
				my $pos = $stop;
				print OUT "$chr\t$pos\t$pos\t.\t0\t$strand\n";
			}
		}else{
			if ($mm9_strand eq "+"){
				#start of annotation in mm9 (+) --> 1st position mm9 ($start)
				#start of annotation in mm10(-) --> 2nd position mm10 ($stop)
				my $pos = $stop;
				print OUT "$chr\t$pos\t$pos\t.\t0\t$strand\n";
			}else{
				#start of annotation in mm9 (-) --> 2nd position mm9 ($stop)
				#start of annotation in mm10(+) --> 1st position mm10 ($start)
				
				my $pos = $start;
				print OUT "$chr\t$pos\t$pos\t.\t0\t$strand\n";
			}
		}#end else
	}#end while (<IN>)
close(IN);

close(OUT);
exit;