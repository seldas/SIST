####################
# SIST: Separate Accugenomics Spike-In reads
# Version 1.01
# 
# #### Updates: ####
# #### End of Updates ###
#
# Contact leihong.wu@fda.hhs.gov for further assistance.
#
####################

use Getopt::Long;
use strict;

########################################################
# USAGE
my $USAGE =<<USAGE;
	Usage:
	  SIST.pl [-fastq_1 <fastq.gz> -fastq_2 <fastq.gz> -O <output_prefix> ] [options] [-help] [-eval]
	
	where:
	-fastq_1  paired-end reads 1 (required)
	-fastq_2  paired-end reads 2 (required)
	-O	Output header (required)
	
	options:
	  -type [match, all, gzip] 
	    match: only extract the spike-in reads 
	    all:   generate fastq file for spike-in reads and remain origin reads (default)
	    gzip:  generate gz file instead of fastq file
	  -t [num]: threads used in BWA mem (default: 8)
	  -bwa_b [num]: mismatch penalty used in BWA mem (default: 6)
	  -bwa_o [num,num]: indel penalty used in BWA mem (default: [20,20])
	  -min_length [num]: minimum match length for a read (default: 70)
	  -clean [true, false]: remove intermediate files (default: true)
	  -filter [true, false]: using human reference genome to improve the result quality (default: true);
	  
	-help:  Prints out this helpful message
	-eval:  Evaluation mode (for research use only, need picard)
	
USAGE
#
######################################################

my $threads = 4;
my $type = 'all'; # in default only do match
my $BWA_B=6; # penalty for mismatch. default = 4 in BWA MEM for short reads
my $BWA_O='20,20';	# penalty for indels. default = 6, high value because indels are highly not expected in spike-in sequence.
my $map_score = 70; # at least mapped 100 bps in consecutive to the spike-in Reference. 
#my $cov = 0; # used for filtering spike-in reads. 
my $clean = 'true';
my $slow = 'true';

GetOptions (   
             "fastq_1=s"    => \my $SAMPLE_1,
			 "fastq_2=s"  	=> \my $SAMPLE_2,
			 "O=s"			=> \my $hq_header,
			 "type=s"		=> \$type, 
			 "t=s"			=> \$threads,
			 "bwa_b=s"		=> \$BWA_B,
			 "bwa_o=s"		=> \$BWA_O,
			 "min_length=s"	=> \$map_score,
			 #"cov=s"		=> \$cov,
			 "clean=s"		=> \$clean,
			 "filter=s"		=> \$slow,
			 q(help)		=> \my $help,
			 q(test_pos)		=> \my $test_pos, # used for development
			 q(test_neg)		=> \my $test_neg, # used for development
			 q(eval)		=> \my $eval
			 );            

if ($help) {
    print "$USAGE\n";
    exit 0;
}

# Test dataset 
if ($test_pos) {
	# test Spike-in Sample
	$SAMPLE_1='data/test_pos_2_1.fastq.gz';
	$SAMPLE_2='data/test_pos_2_2.fastq.gz';
	
	$hq_header = $SAMPLE_1; 
	$hq_header =~ s/_1\.fastq\.gz//g;
	$hq_header =~ s/.*\///g;
	$hq_header = $hq_header;
}
if ($test_neg) {		
	# test 'Negative' Sample
	$SAMPLE_1='data/test_neg_2_1.fastq.gz'; 
	$SAMPLE_2='data/test_neg_2_2.fastq.gz';

	$hq_header = $SAMPLE_1; 
	$hq_header =~ s/_1\.fastq\.gz//g;
	$hq_header =~ s/.*\///g;
	$hq_header = $hq_header;
}

## Check input availability.

if ($type ne 'match' && $type ne 'all' && $type ne 'gzip' ){
	print('Wrong type! should be one in [match, all, gzip] '."\n");
	exit 0;
}

if (! -f 'Refs/picard.jar' && $eval){
	print('Missing Picard. create a link of picard  Refs/picard.jar'."\n");
	exit 0;
}

if (! -f 'Refs/genome.fa' && $slow ne 'false'){
	print('Missing human reference genome file. create a link of human reference genome (e.g., hg38) as Refs/genome.fa'."\n");
	exit 0;
}

if (! -f 'Refs/genome.dict' && $slow ne 'false' && $eval){
	print('Missing human reference genome dictionary. create a link of human reference genome dictionary (by picard) as Refs/genome.dict'."\n");
	exit 0;
}

if (! -f 'Refs/genome.fa.bwt' && $slow ne 'false'){
	print('Missing human reference genome index. create a link of human reference genome dictionary (by BWA index) such as Refs/genome.fa.bwt, etc.'."\n");
	exit 0;
}
	
if (! $SAMPLE_1 ||! $SAMPLE_2 ||! -f $SAMPLE_1 ||! -f $SAMPLE_2 ||! $hq_header){	
	print "Missing essential input files!\n";
	if (! $SAMPLE_1) {
		print "  FASTQ_1 file is not defined!\n" ;
	}elsif (! -f $SAMPLE_1){
		print "  FASTQ_1 file does not exist. Check if the input file name: ".$SAMPLE_1."\n";
	}
	
	if (! $SAMPLE_2) {
		print "  FASTQ_2 file is not defined!\n" ;
	}elsif (! -f $SAMPLE_2){
		print "  FASTQ_2 file does not exist. Check if the input file name: ".$SAMPLE_2."\n";
	}

	print "  Output file (prefix) is not defined.\n" if (! $hq_header);
	print "Use -help to get more information.\n";
	exit 0;
}

if ($SAMPLE_1 !~ /\.fa*s*t*q$/ && $SAMPLE_1 !~ /\.fa*s*t*q\.gz$/){
	print "Not supported Fastq file format for Fastq_1! Must ended with .fastq[.gz] or .fq[.gz] \n";
	exit 0;
}
if ($SAMPLE_2 !~ /\.fa*s*t*q$/ && $SAMPLE_2 !~ /\.fa*s*t*q\.gz$/){
	print "Not supported Fastq file format for Fastq_2! Must ended with .fastq[.gz] or .fq[.gz] \n";
	exit 0;
}
###########CHECK END##########


# Create output files
my $hq_file_pref = $hq_header.".prefil";
my $hq_file = $hq_header.".txt";
my $sam_file = $hq_header.".sam";

# those file only generated in all or gzip mode.
my $spike_fastq_1 = $hq_header."_spikein_1.fastq";
my $spike_fastq_2 = $hq_header."_spikein_2.fastq";
my $origin_fastq_1 = $hq_header."_origin_1.fastq";
my $origin_fastq_2 = $hq_header."_origin_2.fastq";

# if output file already exists; overwrite it.
print "Warning: ".$hq_file." already exists; will overwrite! \n" if (-f $hq_file);

#################Default Parameter Settings##########
# Human reference genome:
my $REF_homo = 'Refs/genome.fa';

# Accugenomics spike-in Reference
my $REF='Refs/Accugenomics_Spikein.fasta';

# Accugenomics spike-in mutations (vcf file)
my $mut_REF='Refs/muts_filtered_sp.vcf';

my %muts=();
open(FH,$mut_REF);
while(<FH>){
	chomp;
	my @array=split("\t");
	$muts{$array[0].":".$array[1]}=1;
}
close FH;
######################################################

open(OFH_SAM_SPK,'>'.$sam_file);

###### STEP1: detect Spike-in Reads ######
my $stop_sign = 0; # early stop for testing larger files
my %high_reads=();
my %doubt_sam=();
my %rep=();

my %reads_count_total=();
my %reads_count_bwa=();
my %reads_count_spk=();

###### spk statistics ######
my %spk_count=();
$spk_count{'>=4'}=0;
$spk_count{'3'}=0;
$spk_count{'2'}=0;
$spk_count{'1'}=0;
			
#print ("Finding Spike-in Reads ...\n");
	
# Find all mapped reads (candidate spike-in reads) from original read set.
open(FH,'bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF.' '.$SAMPLE_1.' '.$SAMPLE_2.' 2>>tmp/bwa_run.log|');
	
open(OFH_h,'>',$hq_file_pref);
print OFH_h "Status\tTotal\tspike_in\tAbsent_Spikein\tOther_Muts\tReads\n";
	
while(<FH>){
	if (/^@/){
		print OFH_SAM_SPK $_;		
	}else{
		my $read_stat = 0; # determine read is spike-in or original.
		my $stop_sign++;
					
		my $line=$_;
		chomp($line);
		my @array=split("\t",$line);
		$reads_count_total{$array[0]}=1;
		if ($array[2] ne '*' && $line =~ /\tMD:Z:[\dACGT]+\t/){	# sequence has to be mapped
			
		##### Pre-Step : reads preparation #####
			$reads_count_bwa{$array[0]}=1;
			my $cigar = $array[5];
			my $sum_length = length($array[9]);
			
			##             INDEL detection  																###
			# In current version, indel are not allowed to happen. They might be in Doubt reads.          
			my $ins_count = 0;
			my $ins_prev = 0;
			if ($cigar =~ /(\d+)M(\d+)I/){
				$ins_prev = $1;
				$ins_count = $2; 
			}
			my $del_count = 0;
			my $del_prev = 0;
			if ($cigar =~ /(\d+)M(\d+)D/){
				$del_prev = $1; 
				$del_count = $2; 
			}
							
			my $flag_indel = 0;
			$flag_indel = 1 if ($ins_count > 0 || $del_count>0);
			###################################################################################################
									
				
			##            FLAG clip : if two-sides softclips												####
			my $left_clip=0;
			my $left_clip_hard=0;
			my $right_clip=0;
			my $right_clip_hard=0;
			$left_clip = $1 if ($cigar =~ /(\d+)S\d+M/); # 
			$right_clip = $1 if ($cigar =~ /\d+M(\d+)S/); # 
			$left_clip_hard = $1 if ($cigar =~ /(\d+)H\d+M/); # 
			$right_clip_hard = $1 if ($cigar =~ /\d+M(\d+)H/); # 
			
			my $flag_clip = 0;
			
			$flag_clip = 1 if ($right_clip >0 || $left_clip >0); # single-side softclip reads 
			$flag_clip = 2 if ($right_clip_hard >0 || $left_clip_hard >0); # single-side hardclip reads 			
			$flag_clip = 5 if ($left_clip>0 && $right_clip>0); # two-side softclip (one short) reads 
			$flag_clip = 6 if ($left_clip_hard>0 && $right_clip_hard>0); # two-side hardclip (one short) reads 
			$flag_clip = 11 if ($left_clip >20 && $right_clip >20); # two-side long softclip reads 
			$flag_clip = 12 if ($left_clip_hard >20 && $right_clip_hard >20); # two-side long hardclip reads
			####################################################################################################
			
			
			## 			  FLAG match length: 																###
			my $mapped_seq = $sum_length-$left_clip-$right_clip;
			
			my $flag_length = 0;	
			$flag_length = ($mapped_seq-$map_score)/10;
			###################################################################################################
			
			
			## 			   Spike-in mutation site															###
			my $muts_included=0;
			foreach($array[3]..($array[3]+$mapped_seq-1)){
				my $id_tmp = $array[2].':'.$_;
				if (exists $muts{$id_tmp}){
					$muts_included ++ ;
				}
			}
			###################################################################################################

			
			####### find spike-in counts		
			if ($flag_indel ==0 && $muts_included >0 ){ #No indels
				
				# Initialize parameters
				my $remain_spikein = 0;
				my $count_ref=0;
				my $muts_rep=$ins_count; #
				
				my $md_str = $1 if $line =~ (/\tMD:Z:([\dACGT]+)\t/);
				my @md_array = split(/[ACGT]/,$md_str);
				my $pos_start = $array[3];
				pop(@md_array);
				foreach(@md_array){
					my $pos_curr = $pos_start+$_;
					$pos_start = $pos_start+$_+1;
					my $sig=0;
					if (exists $muts{$array[2].':'.$pos_curr}){
						$count_ref++;
					}
				}
				$remain_spikein = $muts_included-$count_ref;
				$muts_rep = $muts_rep+scalar(@md_array)-$count_ref;
		
		##### END OF PRE - STEP #####
		
		##### More than 4 Spike-in mutation detected: M>=4 #####
				my $muts_thres = int($mapped_seq/50)+1;	
				if ($remain_spikein>=4 && $muts_rep <$remain_spikein){
					$spk_count{'>=4'} ++ ;
					if ($muts_rep < $muts_thres || $flag_clip < 5 || $flag_length>0){
						if (exists $high_reads{$array[0]}){								
							print OFH_h $high_reads{$array[0]};
							print OFH_h "M4\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
							delete $high_reads{$array[0]};		
						}else{
							$high_reads{$array[0]}="M4\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
							
							$reads_count_spk{$array[0]}=1;
						}
					}
					# most M4 reads doesn't go to human reference comparison.
					$doubt_sam{$line."\n"}=1 if ($muts_rep > $muts_thres);
					
				}elsif ($flag_clip <10){ #No two side (long) soft-clip
					$spk_count{$remain_spikein}++;	
					#### 3 Spike-in Mutation detected: M=3 ####
					if ($remain_spikein ==3){
						# C1: Perfect match 
						if ($count_ref <=1 && $muts_rep+$count_ref < $muts_thres && $flag_length>0){ # (worst case) e.g., 4-3-1-1; 3-3-0-2;
							if (exists $high_reads{$array[0]}){								
								print OFH_h $high_reads{$array[0]};
								print OFH_h "M3\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								delete $high_reads{$array[0]};
								
							}else{
								$high_reads{$array[0]}="M3\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								
								$reads_count_spk{$array[0]}=1;
							}
							
							$doubt_sam{$line."\n"}=1;
						}
					#### 2 spike-in Mutation detected: M=2 ###
					}elsif ($remain_spikein ==2){
						my $sig_M2 = 0;
						$sig_M2 = 1 if ($count_ref ==0 && $muts_rep <=1 && $flag_length>0); # (worst case) e.g., 2-2-0-1
						$sig_M2 = 0 if ($muts_rep ==1 && $flag_clip >=5); # if 2-2-0-1, two-side clip not accepted.
						if ($sig_M2 == 1){
							if (exists $high_reads{$array[0]}){								
								print OFH_h $high_reads{$array[0]};
								print OFH_h "M2\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								delete $high_reads{$array[0]};
								
							}else{
								$high_reads{$array[0]}="M2\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								$reads_count_spk{$array[0]}=1;
							}
							
							$doubt_sam{$line."\n"}=1;
						}
					}
				}
			}
		}
	}
}

#### Unpaired Reads ####
foreach (keys %high_reads){
	print OFH_h $high_reads{$_};
}
########################

close OFH_h;
close FH;

# write unpaired reads to fastq
open(OFH_txt,'>'.$hq_header.".unpair.txt");
open(OFH_fastq,'>'.$hq_header.".unpair.fastq");
foreach (keys %high_reads){
	my $line =$high_reads{$_};
	my @array = split("\t",$line);
	
	print OFH_txt $line;
	
	print OFH_fastq '@'.$array[5]."\n";
	print OFH_fastq $array[14]."\n";
	print OFH_fastq '+'."\n";
	print OFH_fastq $array[15]."\n";	
}
close OFH_txt;
close OFH_fastq;

## Generate paired SAM FILE. ####
foreach (keys %doubt_sam){
	print OFH_SAM_SPK $_;
}
close OFH_SAM_SPK;


##### SLOW MODE: use human reference genome on non-perfect match (match & fair) reads
my $slow_fq_file = $hq_header.".filter.1.fastq";
my $slow_fq_file_2 = $hq_header.".filter.2.fastq";
# my $slow_output = $hq_header.".Href.txt";
my %suspect_reads=();
my %spike_in_reads=();
if ($slow ne 'false'){
	system('java -jar Refs/picard.jar SamToFastq I='.$sam_file." FASTQ=".$slow_fq_file." SECOND_END_FASTQ=".$slow_fq_file_2.
			" VALIDATION_STRINGENCY=LENIENT 2>>tmp/samtofastq.log");
	if (-s $slow_fq_file > 0){
		open(FH,'bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF_homo.' '.$slow_fq_file.' '.$slow_fq_file_2.' 2>>tmp/bwa_run_2.log|');
		# open(OFH,'>'.$slow_output);
		
		while(<FH>){
			if ($_ !~ /^@/){
				my $line=$_;
				chomp($line);
				my @array=split("\t");
				if ($line =~ /MD:Z:([\dATCG]+)\t/){
					my $sum_length_ref= length($array[9]);
					my $md_string = $1;
					my $mut_ref = scalar(split(/[ACGT]/,$md_string));
					
					if ($md_string !~ /[ATCG]0[ATCG]/){
						my $cigar_ref = $array[5];
						my $left_clip=0;
						my $right_clip=0;
						$left_clip = $1 if ($cigar_ref =~ /(\d+)S\d+M/); # 
						$right_clip = $1 if ($cigar_ref =~ /\d+M(\d+)S/); # 
						my $match_ref = $sum_length_ref-$left_clip-$right_clip;
						
						
						my @array_spk = split("\t",$doubt_sam{$array[0]});
						my $sum_length_spk = length($array_spk[9]); 
						
						my $cigar_spk = $array_spk[5];
						my $left_clip=0;
						my $right_clip=0;
						$left_clip = $1 if ($cigar_spk =~ /(\d+)S\d+M/); # 
						$right_clip = $1 if ($cigar_spk =~ /\d+M(\d+)S/); # 
						my $match_spk = $sum_length_spk-$left_clip-$right_clip;
								
						my $md_spk = $1 if ($doubt_sam{$array[0]} =~ /MD:Z:([\dATCG]+)\t/);
						my $mut_spk = scalar(split(/[ACGT]/,$md_spk));
							
							
						if ($match_spk < $match_ref -10 || $md_spk =~ /[ATCG]0[ATCG]/ || $mut_spk > $mut_ref || ($mut_ref ==0 && $cigar_ref =~ /^\d+M$/) ){
							# print OFH "Suspect_reads\t".$line."\n";
							# print OFH "Spike in Map\t".$doubt_sam{$array[0]}."\n";
							$suspect_reads{$array[0]}=1;
						}
					}
				}
			}
		}
		close FH;
	}
	open(FH,'bwa mem -v 0 -SP -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF_homo.' '.$hq_header.'.unpair.fastq 2>>tmp/bwa_run_3.log|');
	while(<FH>){
		if ($_ !~ /^@/){
			my $line=$_;
			chomp($line);
			my @array=split("\t");
			if ($line =~ /MD:Z:([\dATCG]+)\t/){
				my $sum_length_ref = length($array[9]); 
				my $md_string = $1;
				my $mut_ref = scalar(split(/[ACGT]/,$md_string));
				
				if ($md_string !~ /[ATCG]0[ATCG]/){
				
					my $cigar_ref = $array[5];
					my $left_clip=0;
					my $right_clip=0;
					$left_clip = $1 if ($cigar_ref =~ /(\d+)S\d+M/); # 
					$right_clip = $1 if ($cigar_ref =~ /\d+M(\d+)S/); # 
					my $match_ref = $sum_length_ref-$left_clip-$right_clip;
					
					
					my @array_spk = split("\t",$doubt_sam{$array[0]});
					my $sum_length_spk = length($array_spk[9]); 
					
					my $cigar_spk = $array_spk[5];
					my $left_clip=0;
					my $right_clip=0;
					$left_clip = $1 if ($cigar_spk =~ /(\d+)S\d+M/); # 
					$right_clip = $1 if ($cigar_spk =~ /\d+M(\d+)S/); # 
					my $match_spk = $sum_length_spk-$left_clip-$right_clip;
							
					my $md_spk = $1 if ($doubt_sam{$array[0]} =~ /MD:Z:([\dATCG]+)\t/);
					my $mut_spk = scalar(split(/[ACGT]/,$md_spk));
						
						
					if ($match_spk < $match_ref -10 || $md_spk =~ /[ATCG]0[ATCG]/ || $mut_spk > $mut_ref || ( $mut_ref ==0 && $md_string =~ /^\d+M$/ )){
						# print OFH "Suspect_reads\t".$line."\n";
						# print OFH "Spike in Map\t".$doubt_sam{$array[0]}."\n";
						$suspect_reads{$array[0]}=1;
					}
				}
			}
		}
	}	
	close FH;
	# close OFH;
	
	
	open(OFH,'>'.$hq_file);
	print OFH "Total Statistics:\n";
	print OFH "1. Total Output Reads:\t".scalar(keys %reads_count_total)."\n";
	print OFH "2. Total Mapped Reads:\t".scalar(keys %reads_count_bwa)."\n";
	print OFH "3. Total Spike-in Reads (Before Filtering):\t".scalar(keys %reads_count_spk)."\n";
	print OFH "4. Spike-in position (>=4):\t".$spk_count{'>=4'}."\n";
	print OFH "5. Spike-in position (3):\t".$spk_count{'3'}."\n";
	print OFH "6. Spike-in position (2):\t".$spk_count{'2'}."\n";
	print OFH "7. Spike-in position (1):\t".$spk_count{'1'}."\n";

	open(FH,$hq_file_pref);
	while(<FH>){
		my $line = $_;
		my @array=split("\t",$line);
		if (!exists $suspect_reads{$array[5]}){
			print OFH $line;
			$spike_in_reads{$array[5]}=1;
		}
	}
	close FH;	
	close OFH;
		
	
	##### CLEAN UP #####
	if ($clean eq 'true'){
		system('rm -f '.$slow_fq_file);
		system('rm -f '.$slow_fq_file_2);	
		system('rm -f '.$hq_header.'.unpair.fastq');
		system('rm -f '.$hq_header.'.unpair.txt');
		system('rm -f '.$sam_file);
		system('rm -f '.$hq_file_pref);
		# system('rm -f '.$hq_file.'.remained');
		# system('rm -f '.$hq_file.'.filtered');
		# system('rm -f '.$slow_output);
	}
}


##### generate fastq files for spike-in and origin reads #####
if ($type eq 'all' || $type eq 'gzip'){
	print "Warning: ".$spike_fastq_1." already exists; will overwrite! \n" if (-f $spike_fastq_1);
	print "Warning: ".$spike_fastq_2." already exists; will overwrite! \n" if (-f $spike_fastq_2);
	print "Warning: ".$origin_fastq_1." already exists; will overwrite! \n" if (-f $origin_fastq_1);
	print "Warning: ".$origin_fastq_2." already exists; will overwrite! \n" if (-f $origin_fastq_2);
	
	# separating fastq file 1
	print ('Extracting Reads... (Pair - 1)'."\n");
	open(FH,'gzip -cd '.$SAMPLE_1.'|');
	open(OFH_spike,'>'.$spike_fastq_1);
	open(OFH_origin,'>'.$origin_fastq_1);
	my $signal = 0;
	my $count_pair_1=0;
	while(<FH>){
		my $line = $_;
		if (/^@(\S+)[\s\/]+(\S+)/){
			my $array=$1;
			if (exists $spike_in_reads{$array}){
				$count_pair_1++;
				$signal = 1;
			}else{
				$signal = 0;
			}
		} 
		print OFH_spike $line if ($signal ==1);
		print OFH_origin $line if ($signal ==0);
	}
	close FH;
	close OFH_spike;
	close OFH_origin;
	
	# separating fastq file 2
	print ('Extracting Reads... (Pair - 2)'."\n");
	open(FH,'gzip -cd '.$SAMPLE_2.'|');
	open(OFH_spike,'>'.$spike_fastq_2);
	open(OFH_origin,'>'.$origin_fastq_2);
	my $signal = 0;
	# my $count_pair_2=0;
	while(<FH>){
		my $line = $_;
		if (/^@(\S+)[\s\/]+(\S+)/){
			my $array=$1;
			if (exists $spike_in_reads{$array}){
				# $count_pair_2++;
				$signal = 1;
			}else{
				$signal = 0;
			}
		} 
		print OFH_spike $line if ($signal ==1);
		print OFH_origin $line if ($signal ==0);
	}
	close FH;
	close OFH_spike;
	close OFH_origin;
	
	# End of the processing.
	print ('Done! In total '.$count_pair_1." (pairs of) Spike-in reads have been separated.\n");
}

######   Generating Gzip files ######
if ($type eq 'gzip'){
	print ('Gzip output files.. (may take long time...)'."\n");
	system('gzip -1 '.$origin_fastq_1);
	system('gzip -1 '.$origin_fastq_2);
	system('gzip '.$spike_fastq_1);
	system('gzip '.$spike_fastq_2);
	print ('Program finished normally.'."\n");
}

###### evaluation on the separation file ######
if ($eval){
	#measure insert size in separation reads	
	system('java -jar Refs/picard.jar CollectInsertSizeMetrics I='.$sam_file.' O=insert_spikein.txt H=insert_spikein.pdf M=0.5');
}

###### CLEAN UP ######
if ($clean eq 'true'){
	#system('rm -f tmp/bwa_run.log'); # remove BWA log file 
	#system('rm -f tmp/bwa_run_2.log'); # remove BWA log file 
	#system('rm -f tmp/bwa_run_3.log'); # remove BWA log file 
	#system('rm -f tmp/samtofastq.log'); # remove BWA log file 
}
