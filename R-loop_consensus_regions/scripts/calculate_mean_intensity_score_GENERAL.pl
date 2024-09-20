#!/usr/bin/perl



&calculate_mean_intensity_score;

sub calculate_mean_intensity_score
{
	my $file_in = $ARGV[0];
	open(INFO, $file_in);

	my $file_out = $file_in."_after_computing_average.bed";
	open(OUT, ">$file_out");

	while(<INFO>)
	{
		chomp;
		my ($chr, $start, $end, $scores) = split(/\t/,$_,4);
		my @cols = split(/\t/,$scores);
		my $total = 0;
		my $count = 0;
		for(my $i =0; $i<@cols; $i++)
		{
			$total = $total + $cols[$i];
			$count++;
		}
		my $average = $total/$count;
		print OUT "$chr\t$start\t$end\t$average\n";
	}
	close(INFO);
	close(OUT);
}
