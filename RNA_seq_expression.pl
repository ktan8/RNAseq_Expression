#!/usr/bin/perl

use strict;
use warnings;
use File::Path qw(make_path);
use Parallel::ForkManager;

# The following script is a script written to perform the RPKM
# based measurement of the expression of the gene
#my $annotation_file = "/home/kartong/Google\ Drive/Code/RPKM_measurement/hg19_refseq_20140619_mod.txt";
#my $output_folder = "/home/kartong/ouput_folder/";
#my $bamfile_list = "/home/kartong/sample_list.txt.csv";

# Input arguments
my $annotation_file = $ARGV[0];
my $output_folder 	= $ARGV[1];
my $bamfile_list 	= $ARGV[2];

my $file_info_list		= $output_folder . "RNAseq_expr_bam_info.txt";
my $combined_rpkm_list	= $output_folder . "RNAseq_expr_combined_rpkm.txt";


my $pm = Parallel::ForkManager->new(8);


# Parse the information I need from the annotation table
my ($data_AOA, $header_hash) = parse_refseq_annotation($annotation_file);
make_path($output_folder);


my @overall_info_AOA;
my %overall_RPKM_HOA;


$pm -> run_on_finish ( # called BEFORE the first call to start()
	sub {
		my ($pid, $exit_code, $ident, $exit_signal, $core_dump, $data_structure_reference) = @_;
		my $rpkm_array_ref 	= $data_structure_reference -> [0];
		my $total_reads 	= $data_structure_reference -> [1];
		my $mapped_reads	= $data_structure_reference -> [2];
		my $read1			= $data_structure_reference -> [3];	
		my $read2			= $data_structure_reference -> [4];
		my $label			= $data_structure_reference -> [5];
		
		
		# Store the information somewhere
		my $overall_read_info = [$label, $total_reads, $mapped_reads, $read1, $read2];
		push(@overall_info_AOA, $overall_read_info);
		$overall_RPKM_HOA{ $label } = $rpkm_array_ref;
	}
);





open(INPUT, $bamfile_list) || die $!;

# Loop through a list of bam files where
# each line is a label, followed by the 
# bamfile path
while(my $line = <INPUT>){
	$line =~ s/\r\n//;
	chomp($line);
	my ($label, $bamfile) = split(/\t/, $line);
	
	$pm -> start() and next;
	
	my $outputfile = $output_folder . $label . ".expression.txt";
	
	open(my $EXPR, ">", $outputfile);

	print $EXPR join("\t","transcript_id", "genesymbol", "readcount", "transcript_length", "RPKM") . "\n";

	my ($total_reads, $mapped_reads, $read1, $read2) = count_mapped_reads($bamfile);
	#my $mapped_reads = 248035594;


	# Loop through all the annotations
	my @RPKM_array;
	foreach(@$data_AOA){
		my $chr				= $_ ->[$header_hash -> {chrom}	];
		my $start			= $_ ->[$header_hash -> {txStart}];
		my $end				= $_ ->[$header_hash -> {txEnd}	];
		my $exonstart_str	= $_ ->[$header_hash -> {exonStarts}];
		my $exonend_str		= $_ ->[$header_hash -> {exonEnds}	];
		my $genesymbol		= $_ ->[$header_hash -> {name2}	];
		my $transcript_id	= $_ ->[$header_hash -> {name}	];
	
	
		my $readcount 	= count_reads_in_region($chr, $start, $end, $bamfile);
		my $genelength 	= calculate_gene_length($exonstart_str, $exonend_str);
	

		# Get RPKM value
		my $RPKM = calculate_rpkm($readcount, $genelength, $mapped_reads);
	
		# Print one list for each sample
		print $EXPR join("\t", $transcript_id, $genesymbol, $readcount, $genelength, $RPKM) . "\n";
		
		push(@RPKM_array, $RPKM);
	}
		
	close($EXPR);
	
	my $result = [\@RPKM_array, $total_reads, $mapped_reads, $read1, $read2, $label];
	
	
	$pm->finish(0, $result);
}

$pm->wait_all_children;

close(INPUT);


#====================================================================

# Print out the general result
# Also print a combined list for all samples




open(FILEINFO, ">", $file_info_list);


print FILEINFO join("\t", "label", "total_reads", "mapped_reads", "read1", "read2") . "\n";
foreach(@overall_info_AOA){
	print FILEINFO join("\t", @{ $_ }) . "\n";
}
close(FILEINFO);


# Print the output to a combined
# RPKM file
open(COMB_RPKM, ">", $combined_rpkm_list);
print COMB_RPKM join("\t", "transcript_id", "gene_name", keys(%overall_RPKM_HOA) ) . "\n";
my $sample_label;
foreach(keys(%overall_RPKM_HOA)){ $sample_label = $_ };


my $transcript_num = scalar( @{ $overall_RPKM_HOA{$sample_label} } );


for(my $i=0; $i< $transcript_num; $i++){
	my @line_rpkm;
	foreach(values(%overall_RPKM_HOA)){
		my $sample_rpkm = $_ -> [$i];
		push(@line_rpkm, $sample_rpkm);
	}

	my $line_ref 		= $data_AOA->[$i];
	my $genesymbol		= $line_ref ->[$header_hash -> {name2}	];
	my $transcript_id	= $line_ref ->[$header_hash -> {name}	];

	print COMB_RPKM join("\t", $transcript_id, $genesymbol, @line_rpkm) . "\n";

}
close(COMB_RPKM);



#====================================================================


sub parse_refseq_annotation{
	my $ref_seq_file = $_[0];

	open(my $REFSEQ, $ref_seq_file) || die $!;
	my $header = <$REFSEQ>;
	$header =~ s/\r\n//;
	chomp($header);

	# Convert header information into array number
	$header =~ s/^\#//;
	my @header_array = split(/\t/, $header);
	my %header_hash;
	
	my $counter = 0;
	foreach(@header_array){
		$header_hash{$_} = $counter;
		$counter++;
	}
	
	# Store all data into AOA
	my @AOA;
	while(my $line = <$REFSEQ>){
		$line =~ s/\r\n//;
		chomp($line);
		
		my @line_array = split(/\t/, $line);
		push(@AOA, \@line_array);
	}
		
	close($REFSEQ);
	
	
	return(\@AOA, \%header_hash);
}


#====================================================================


# Use samtools to count number of reads in the region of interest
sub count_reads_in_region{
	my $chr 		= $_[0];
	my $start		= $_[1];
	my $end			= $_[2];
	my $bam_file	= $_[3];


	my $region_str = $chr . ":" . $start . "-" . $end;
	
	#print "samtools view $bam_file $region_str | wc -l" . "\n";
	
	my $read_count = `samtools view $bam_file $region_str | wc -l`;
	chomp($read_count);
	
	
	return ($read_count);
}




#====================================================================


sub calculate_gene_length{
	my $exonstart_str	= $_[0];
	my $exonend_str		= $_[1];

	my @start_array	= split(",", $exonstart_str);
	my @end_array 	= split(",", $exonend_str);
	
	my $genelength =0;
	
	foreach(my $i=0; $i<scalar(@start_array); $i++){
		my $start 	= $start_array[$i];
		my $end		= $end_array[$i];
		my $length 	= $end - $start;
		$genelength += $length;
	}
	
	return $genelength;
}




#====================================================================


sub count_mapped_reads{
	my $bamfile = $_[0];
	
	my @flagstat_result = split("\n", `samtools flagstat $bamfile`);

	my $total_reads		= ( $flagstat_result[0] =~ m/^([0-9]+)/ )[0];
	my $mapped_reads	= ( $flagstat_result[2] =~ m/^([0-9]+)/ )[0];
	my $read1			= ( $flagstat_result[4] =~ m/^([0-9]+)/ )[0];
	my $read2			= ( $flagstat_result[5] =~ m/^([0-9]+)/ )[0];
	
	#284052914 + 0 in total (QC-passed reads + QC-failed reads)
	#0 + 0 duplicates
	#248035594 + 0 mapped (87.32%:-nan%)
	#284052914 + 0 paired in sequencing
	#142060388 + 0 read1
	#141992526 + 0 read2
	#223121430 + 0 properly paired (78.55%:-nan%)
	#234505274 + 0 with itself and mate mapped
	#13530320 + 0 singletons (4.76%:-nan%)
	#5898040 + 0 with mate mapped to a different chr
	#5898022 + 0 with mate mapped to a different chr (mapQ>=5)

	
	return ($total_reads, $mapped_reads, $read1, $read2);
}


#====================================================================


sub calculate_rpkm{
	my $genecounts	= $_[0];
	my $genelength	= $_[1];
	my $library_size= $_[2];

	my $RPKM = $genecounts / ($genelength * $library_size) * 10 ** 9;

	return $RPKM
}



