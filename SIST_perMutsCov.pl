$bed_file_ref = 'Refs/muts_filtered_ref_hg37.bed';
$bed_file_spk = 'Refs/muts_filtered_sp.bed';

$origin_folder = 'Thermo_cfDNA';
$spike_folder = 'result_Thermo_cfDNA';

opendir(DH, $origin_folder);
@origin_folder_sub = readdir(DH);
close DH;

system('mkdir -p Cov.'.$spike_folder);
foreach(@origin_folder_sub){
	if ($_ =~ /^Nala/){  # only fDNA + ACG are analyzed.
		$sub_folder = $_;
		opendir(DH_sub, $origin_folder.'/'.$sub_folder);
		@raw_file = readdir(DH_sub);
		close DH_sub;
		
		foreach(@raw_file){
			if (/(Tag.*)\.bam$/){
				$file_header = $1;
				$test_bam_ref = $spike_folder.'/'.$sub_folder.'.'.$file_header.'.b6.l70_origin_1.bam';
				$test_sam_spk = $spike_folder.'/'.$sub_folder.'.'.$file_header.'.b6.l70.sam';
				$output_file  = 'Cov.'.$spike_folder.'/'.$sub_folder.'.'.$file_header.'.cov.txt';
				
				getCov($test_bam_ref,$test_sam_spk, $output_file);
				
			}
		}
	}
}

sub getCov($test_bam_ref,$test_sam_spk, $output_file){

 	## spike-in bam ##
	$test_sam_spk_name = $test_sam_spk;
	$test_sam_spk_name =~ s/.+\///g;
	$tmp_bam = 'tmp/'.$test_sam_spk_name;
	$tmp_bam =~ s/\.sam/\.bam/g;

	%cover_spk=();

	system('samtools view -b '.$test_sam_spk.' |samtools sort -o '.$tmp_bam.' ');
	system('samtools index '.$tmp_bam);

	%spk_read_sequence=();
	open(FH_spk,'samtools view '.$tmp_bam.'|');
	while(<FH_spk>){
		chomp;
		@array = split("\t");
		$spk_read_sequence{$array[0]}=1;
	}
	close FH_spk;

	$i=0;
	open(FH,$bed_file_spk);
	while(<FH>){
		chomp;
		$line = $_;
		$cover_spk{$i}=0;	
		open(FH_spk,'samtools view '.$tmp_bam.' '.$line.'|');
		%tmp_readname=();
		while(<FH_spk>){
			chomp;
			@array=split("\t");
			if (!exists $tmp_readname{$array[0]}){
				$cover_spk{$i} ++ ;
				$tmp_readname{$array[0]}=1;
			}
		}
		$cover_spk{$i}=$line."\t".$cover_spk{$i};
		close FH_spk;
		$i++;
	}

	## Mutation positions ##

	#unless (-e $test_bam_ref.'.bai'){
		system('samtools index '.$test_bam_ref);
	#}
	%cover_ref=();
	$i=0;
	open(FH,$bed_file_ref);
	while(<FH>){
		$line = $_;
		chomp($line);
		$cover_ref{$i}=0;
		open(FH_ref,'samtools view '.$test_bam_ref.' '.$line.'|');
		%tmp_readname=();
		while(<FH_ref>){
			chomp;
			@array = split("\t");
			if (!exists $tmp_readname{$array[0]}){
				$cover_ref{$i}++ if (!exists $spk_read_sequence{$array[0]});
				$tmp_readname{$array[0]}=1;
			}
		}
		$cover_ref{$i} = $line."\t".$cover_ref{$i};
		close FH_ref;
		$i++;
	}
	close FH;

	open(OFH,'>'.$output_file);
	print OFH "Input:\n";
	print OFH $test_bam_ref."\n";
	print OFH $test_sam_spk."\n\n";
	print OFH "Reference region\tReference Count\tSpike-in Region\tSpike-in Count\n";
	foreach(sort keys %cover_spk){
		@array_spk = split("\t",$cover_spk{$_});
		print OFH $cover_ref{$_}."\t".$cover_spk{$_}."\n" if ($array_spk[1]> 0);
	}
	close OFH;

	print "Finshed! ".$output_file."\n";
} 


 

