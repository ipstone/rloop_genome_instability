#!/usr/bin/perl


&loop_through_files;


sub loop_through_files
{
	my $file_list = "list_of_merged_strand_files.txt";
	open(LIST, $file_list);

	while(<LIST>)
	{
		chomp;
		my $file = $_.".bed";
		&calculate_max_intensity_score($file);
	}
	close(LIST);
}

sub calculate_max_intensity_score
{
	my $file_in = shift;
	open(INFO, $file_in);

	my $file_out = $file_in."_after_getting_max_across_strands.txt";
	open(OUT, ">$file_out");

	while(<INFO>)
	{
		chomp;
		my ($chr, $start, $end, $strand1, $strand2) = split(/\t/,$_,5);
		my $max = $strand1;
		if($strand1 > $strand2)
		{
			$max = $strand1;
		}
		if($strand2 > $strand1)
                {
                        $max = $strand2;
                }

		print OUT "$chr\t$start\t$end\t$max\n";
	}
	close(INFO);
	close(OUT);
}
