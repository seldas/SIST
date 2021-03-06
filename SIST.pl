####################
# SIST: Separate Accugenomics Spike-In reads
# Version 1.30
# 
# #### Updates: ####
# 08/16/2018: using git-hub for version control (update information depreciated)
# 08/02/2018: provide an option to use customized spike-in reference. (Ver. 1.10 -> 1.20)
# 06/09/2018: Support BAM output (to include additional read information)
# 01/04/2018: Support BAM input (Ver. 1.01 -> 1.10)
# 12/02/2017: Fix a Bug about read head detection in FASTQ file generation step. (ver. 1.00 -> 1.01)
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
      SIST.pl [-fastq_1 <fastq.gz> -fastq_2 <fastq.gz>] [-bam <bam>] -O <output_prefix> [options]
    
    Required inputs:
      -fastq_1  paired-end reads 1 (required if no bam input)
      -fastq_2  paired-end reads 2 (required if paired-end reads)
    
      -bam  input file as bam format (required if no fastq input)
    
      -O    Output file header 
    
    Common options:
      -type [match, all, gzip] 
        match: only extract the spike-in reads 
        all:   generate fastq file for spike-in reads and remain origin reads (default)
        gzip:  generate gz file instead of fastq file
      -keepBam : when the input is bam file and keepBam is specified, output would keep the sam format and not change to fastq; 
      -filter [true, false]: using human reference genome to improve the result quality (default: true);
      -ref [prefix]: prefix of spike-in reference, (defalt: "Refs/Accugenomics_Spikein");
      
    Other options:
      -t [num]: threads used in BWA mem (default: 8);
      -bwa_b [num]: mismatch penalty used in BWA mem (default: 6);
      -bwa_o [num,num]: indel penalty used in BWA mem (default: [20,20]);
      -min_length [num]: minimum match length for a read (default: 70);
      -umi [true, false]: flag for umi reads (default: false);
      -clean [true, false]: remove intermediate files (default: true);
      
    -help:  Prints out this helpful message
    
USAGE
#
######################################################

my $SAMPLE_1 = '';
my $SAMPLE_2 = '';
my $BAM = '';
my $hq_header = '';

my $threads = 4;
my $type = 'all';
my $BWA_B=6; # penalty for mismatch. default = 4 in BWA MEM for short reads
my $BWA_O='20,20';  # penalty for indels. default = 6, high value because indels are highly not expected in spike-in sequence.
my $map_score = 70; # at least mapped 100 bps in consecutive to the spike-in Reference. 
#my $cov = 0; # was used for filtering spike-in reads. 
my $clean = 'true';
my $slow = 'true';
my $umi_tag = 'false';
my $single_mode = 'false';
my $samtools_path = 'samtools';
my $test_mode = 'false';
my $ref_prefix = 'Refs/Accugenomics_Spikein';

#### TESTING MODE ONLY #####
# my $samtools_path = '/storage2/lwu/SEQC2/samtools-1.8/build/bin/samtools'; # specific samtools path for test use.
############################

GetOptions (   
            "fastq_1=s" => \$SAMPLE_1,
            "fastq1=s"      => \$SAMPLE_1,
        "fq_1=s"        => \$SAMPLE_1,
        "fq1=s"     => \$SAMPLE_1,
        "f1=s"      => \$SAMPLE_1,
        "1=s"       => \$SAMPLE_1,
        
        "fastq2=s"      => \$SAMPLE_2,
        "fastq_2=s"     => \$SAMPLE_2,
        "fq_2=s"    => \$SAMPLE_2,
        "fq2=s"     => \$SAMPLE_2,
        "f2=s"      => \$SAMPLE_2,
        "2=s"       => \$SAMPLE_2,
             
        "bam=s"     => \$BAM,
             
        "O=s"       => \$hq_header,
        "o=s"       => \$hq_header,
             
        "ref=s"     => \$ref_prefix,
        "r=s"       => \$ref_prefix,
        "R=s"       => \$ref_prefix,
        "Ref=s"     => \$ref_prefix,
        
        "umi=s"     => \$umi_tag,
             
        "type=s"    => \$type, 
        "t=s"       => \$threads,
        "bwa_b=s"   => \$BWA_B,
        "bwa_o=s"   => \$BWA_O,
        "min_length=s"  => \$map_score,
        #"cov=s"    => \$cov,
        "clean=s"   => \$clean,
        "filter=s"  => \$slow,
        q(help)     => \my $help,
        q(h)        => \my $help,
        q(eval)     => \my $eval,
        q(keepBam)  => \my $keepBam,
        # ONLY USED FOR DEV PURPOSES;
        "test_mode=s"   => \$test_mode
);            

## TEST_MODE                     
($SAMPLE_1, $SAMPLE_2, $BAM, $hq_header) = get_test_data($test_mode) if ($test_mode ne 'false');

## Code Section 1. Pre-Check Begin ##
# c1.1 Help message
if ($help || (!$SAMPLE_1 && !$BAM)) {
    # print "\tSomething Wrong\n";
    print "$USAGE\n";
    exit 0;
}

# c1.2 Check analysis type (-type)
if ($type ne 'match' && $type ne 'all' && $type ne 'gzip' ){
    print('Wrong type! should be one in [match, all, gzip] '."\n");
    exit 0;
}

# c1.3 Check Picard (only works when doing insert size measurement)
#if (! -f 'Refs/picard.jar' && $eval){
#   print('Missing Picard. create a link of picard  Refs/picard.jar'."\n");
#   exit 0;
#}

# c1.4 Check Human Reference Files, only checked when "slow" mode is on 
if (! -f 'Refs/genome.fa' && $slow ne 'false'){
    print('Missing human reference genome file. 
        Create a link of human reference genome (e.g., hg38) as Refs/genome.fa'."\n");
    exit 0;
}

if (! -f 'Refs/genome.dict' && $slow ne 'false' && $eval){
    print('Missing human reference genome dictionary. 
        Create a link of human reference genome dictionary (by picard) as Refs/genome.dict'."\n");
    exit 0;
}

if (! -f 'Refs/genome.fa.bwt' && $slow ne 'false'){
    print('Missing human reference genome index. 
        Create a link of human reference genome dictionary (by BWA index) such as Refs/genome.fa.bwt, etc.'."\n");
    exit 0;
}

# c1.5 Check Input Files Availability (-f1, -f2, -bam)
if ( (!$SAMPLE_2)&& !$BAM){
    $single_mode = 'true'; # single-end reads are used in this study
    print "Detect current analysis uses single-ended reads; 
        if not, please stop and double check whether fastq_2 is correct. \n";
}

if ( (!$BAM||! -f $BAM) && (!$SAMPLE_1||! -f $SAMPLE_1 )){  
    print $BAM ;
    print " Missing input files!\n";
    print " Use -help to get more information.\n";
    exit 0;
}

if ( (!$BAM||! -f $BAM) && ($single_mode eq 'false') && (!$SAMPLE_2||! -f $SAMPLE_2 )){ 
    print " Missing paired input files!\n";
    print " Use -help to get more information.\n";
    exit 0;
}

if ($BAM){
    if ($BAM !~ /\.bam$/ ){
        print "Not supported BAM input file format! Must ended with .bam \n";
        exit 0;
    }
}

if ($SAMPLE_1){
    if ($SAMPLE_1 !~ /\.fa*s*t*q$/ && $SAMPLE_1 !~ /\.fa*s*t*q\.gz$/){
        print "Not supported Fastq file format for Fastq_1! Must ended with .fastq[.gz] or .fq[.gz] \n";
        exit 0;
    }
}

if ($SAMPLE_2){
    if ($SAMPLE_2 !~ /\.fa*s*t*q$/ && $SAMPLE_2 !~ /\.fa*s*t*q\.gz$/){
        print "Not supported Fastq file format for Fastq_2! Must ended with .fastq[.gz] or .fq[.gz] \n";
        exit 0;
    }
}

# c1.6 Check Output Header
if (!$hq_header){   
    print " Missing output folder!\n";
    print " Use -help to get more information.\n";
    exit 0;
}

# c1.7 Create tmp folder if not existed.
my $tmp_dir = 'tmp';
mkdir($tmp_dir) unless(-d $tmp_dir);

## Pre-Check End ##


## Code Section 2. Variable Preparation Begin ##
# c2.1 Create output files
my $hq_file_pref = $hq_header.".prefil";
my $hq_file = $hq_header.".txt";
my $hq_file_suspected = $hq_header."_suspected.txt";
my $sam_file = $hq_header.".sam";

# those file only generated in all or gzip mode.
my $spike_fastq_1 = $hq_header."_spikein_1.fastq";
my $origin_fastq_1 = $hq_header."_origin_1.fastq";
my $suspect_fastq_1 = $hq_header."_suspect_1.fastq";

$spike_fastq_1 = $hq_header."_spikein_1.bam" if ($keepBam);
$origin_fastq_1 = $hq_header."_origin_1.bam" if ($keepBam);
$suspect_fastq_1 = $hq_header."_suspect_1.bam" if ($keepBam);

my $spike_fastq_2 = $hq_header."_spikein_2.fastq" if($SAMPLE_2);
my $origin_fastq_2 = $hq_header."_origin_2.fastq" if($SAMPLE_2);
my $suspect_fastq_2 = $hq_header."_suspect_2.fastq" if($SAMPLE_2);

# if output file already exists; overwrite it.
print "Warning: ".$hq_file." already exists; will overwrite! \n" if (-f $hq_file);

# Human reference genome:
my $REF_homo = 'Refs/genome.fa';

# Accugenomics spike-in Reference
my $REF=$ref_prefix.'.fasta';
my $mut_REF='Refs/muts_filtered_sp.vcf';
if ($ref_prefix ne 'Refs/Accugenomics_Spikein'){
    $mut_REF=prepare_ref($ref_prefix);
}

my %muts=();
open(FH,$mut_REF);
while(<FH>){
    chomp;
    my @array=split("\t");
    $muts{$array[0].":".$array[1]}=1;
}
close FH;

## Variable Preparation End ##


## Code Section 3. Main Program Begin ##
###### c3.1: detect Spike-in Reads ######
## c3.1.1 initialize 
my $stop_sign = 0; # early stop for testing larger files
my %high_reads=();
my %doubt_sam=();
my %doubt_sam_ID=();
my %reliable_sam=();
my %rep=();

my %reads_count_total=();
my %reads_count_bwa=();
my %reads_count_spk=();

# spk statistics
my %spk_count=();
$spk_count{'>=4'}=0;
$spk_count{'3'}=0;
$spk_count{'2'}=0;
$spk_count{'1'}=0;      
#print ("Finding Spike-in Reads ...\n");

# output file (temporary)
open(OFH_h,'>',$hq_file_pref);
print OFH_h "Status\tTotal\tspike_in\tAbsent_Spikein\tOther_Muts\tReads\n";
open(OFH_SAM_SPK,'>'.$sam_file);

# c3.1.2 Find all mapped reads (candidate spike-in reads) from original read set.
if ($SAMPLE_1){
    if ($single_mode eq 'false'){
        # paired-end reads
        open(FH,'bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF.' '.
            $SAMPLE_1.' '.$SAMPLE_2.' 2>tmp/bwa_run.log|');
    }else{
        # single-end reads
        open(FH,'bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF.' '.
            $SAMPLE_1.' 2>tmp/bwa_run.log|');
    }
}elsif($BAM){
    # input is Bam
    open(FH,$samtools_path.' fastq '.$BAM.'|bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF.
            ' - 2>>tmp/bwa_run.log|');
}else{
    exit 0;
}

# read output from BWA alignment
while(<FH>){
	if(/^@/){
		# print header information
		print OFH_SAM_SPK $_;     
	}else{
		# determine read is spike-in or original.
		my $read_stat = 0; 
		my $stop_sign++;

		my $line=$_;
		chomp($line);
		my @array=split("\t",$line);
		$reads_count_total{$array[0]}=1;

		# sequence has to be mapped
		if ($array[2] ne '*' && $line =~ /\tMD:Z:[\d\^ACGT]+\t/){         
			# Pre-Step : reads preparation
			$reads_count_bwa{$array[0]}=1;
			my $cigar = $array[5];
			my $sum_length = length($array[9]);

			# INDEL detection
			my $ins_count = 0;
			my $ins_prev = 0;
			if ($cigar =~ /(\d+)M(\d+)I/){
				$ins_count = $ins_count+$& while($cigar =~ /(\d+)I/g );
				# $ins_prev = $1;
				# $ins_count = $2; 
			}

			my $del_count = 0;
			my $del_prev = 0;
			if ($cigar =~ /(\d+)M(\d+)D/){
				$del_count = $del_count+$& while($cigar =~ /(\d+)D/g );
				# $del_prev = $1; 
				# $del_count = $2; 
			}

			my $flag_indel = 0;
			$flag_indel = 1 if ($ins_count > 0 || $del_count>0);                         

			# FLAG for clipping : if two-sides softclips
			my $left_clip=0;
			my $left_clip_hard=0;
			my $right_clip=0;
			my $right_clip_hard=0;
			$left_clip = $1 if ($cigar =~ /(\d+)S\d+M/); # 
			$right_clip = $1 if ($cigar =~ /\d+M(\d+)S/); # 
			$left_clip_hard = $1 if ($cigar =~ /(\d+)H\d+M/); # 
			$right_clip_hard = $1 if ($cigar =~ /\d+M(\d+)H/); # 
				
			my $flag_clip = 0;
			if ($umi_tag eq 'true'){
				#no clipping criterion for umi reads;
				$flag_clip = 5;
			}else{
				#if not umi reads, soft-clipping criterion still works
				$flag_clip = 1 if ($right_clip >0 || $left_clip >0); # single-side softclip reads 
				$flag_clip = 2 if ($right_clip_hard >0 || $left_clip_hard >0); # single-side hardclip reads             
				$flag_clip = 5 if ($left_clip>0 && $right_clip>0); # two-side softclip (one short) reads 
				$flag_clip = 6 if ($left_clip_hard>0 && $right_clip_hard>0); # two-side hardclip (one short) reads 
				$flag_clip = 11 if ($left_clip >20 && $right_clip >20); # two-side long softclip reads 
				$flag_clip = 12 if ($left_clip_hard >20 && $right_clip_hard >20); # two-side long hardclip reads
			}
				
			## FLAG for maximum match length:
			my $mapped_seq = $sum_length-$left_clip-$right_clip-$ins_count+$del_count;
			my $flag_length = 0; 
			$flag_length = ($mapped_seq-$map_score)/10;

			## counting Spike-in mutation site 
			my $muts_included=0;
			foreach($array[3]..($array[3]+$mapped_seq-1)){
				my $id_tmp = $array[2].':'.$_;
				if (exists $muts{$id_tmp}){
					$muts_included ++ ;
				}
			}
							
			## find spike-in counts      
			#if ($flag_indel ==0 && $muts_included >0 ){ #No indels
			if ($muts_included >0 ){ #allow indels
				# Initialize parameters
				my $remain_spikein = 0;
				my $count_ref=0;
				my $muts_rep=$ins_count; 

				# Read MD string in each read
				my $md_str = $1 if $line =~ (/\tMD:Z:([\^\dACGT]+)\t/);
				my @md_array = split(/[\^ACGT]+/,$md_str);
				my $pos_start = $array[3];
				pop(@md_array);
		
				# Measuring on-site mutations
				foreach(@md_array){
					my $pos_curr = $pos_start+$_;
					$pos_start = $pos_start+$_+1;
					my $sig=0;
					if (exists $muts{$array[2].':'.$pos_curr}){
						$count_ref++;
					}
				}
					
				# Measuring off-site mutations
				$remain_spikein = $muts_included-$count_ref;
				$muts_rep = $muts_rep+scalar(@md_array)-$count_ref; 

				## More than 4 Spike-in mutation detected: M>=4 #####
				my $muts_thres = int($mapped_seq/50)+1; 
			
				if ($remain_spikein>=4 ){ 
					#if on-site mutation (SICs) >4, labeled as high-confidence (M4 group)
					$spk_count{'>=4'} ++ ;
					if ($muts_rep < $muts_thres || $flag_clip < 5 || $flag_length>0){
					# for umi reads, flag_clip did not work.
						if (exists $high_reads{$array[0]}){
							# for reads already paired, output to the final output with priority. 
							print OFH_h $high_reads{$array[0]};
							print OFH_h "M4\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
							delete $high_reads{$array[0]};       
						}else{
							# append the reads into temp pool, waiting for its paired reads.
							$high_reads{$array[0]}="M4\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
							$reads_count_spk{$array[0]}=1;
						}
					}
					#most M4 reads doesn't go through the double check with human reference.
					if ($muts_rep > $muts_thres){
						$doubt_sam{$line."\n"}=1;
						$doubt_sam_ID{$array[0]}=$line."\n";
					}else{
						$reliable_sam{$line."\n"}=1;
					}
					
				}elsif ($flag_clip <10){ 
					#no two side (long) soft-clip, if umi-reads, automatically passed this criterion.
					$spk_count{$remain_spikein}++; 
					if ($remain_spikein == 3){
						## 3 Spike-in Mutation detected: M3
						if ($count_ref <=1 && $muts_rep+$count_ref <= $muts_thres && $flag_length>0){ # (worst case) e.g., 4-3-1-1; 3-3-0-2;
							if (exists $high_reads{$array[0]}){                              
								print OFH_h $high_reads{$array[0]};
								print OFH_h "M3\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								delete $high_reads{$array[0]};
									
							}else{
								$high_reads{$array[0]}="M3\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								$reads_count_spk{$array[0]}=1;
							}
						
							#all M3 reads go through the slow mode check;
							$doubt_sam{$line."\n"}=1;
							$doubt_sam_ID{$array[0]}=$line."\n"; 
						}
					}elsif ($remain_spikein ==2){
						## 2 spike-in Mutation detected: M2
						my $sig_M2 = 0;
						# worst accpeted case: 2-2-0-1
						$sig_M2 = 1 if ($count_ref ==0 && $muts_rep <=1 && $flag_length>0); 
						# if 2-2-0-1, two-side clip not accepted. umi reads automatically passed the rule,
						$sig_M2 = 0 if ($muts_rep ==1 && $flag_clip >5); 
				
						if ($sig_M2 == 1){
							if (exists $high_reads{$array[0]}){                              
								print OFH_h $high_reads{$array[0]};
								print OFH_h "M2\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								delete $high_reads{$array[0]};
									
							}else{
								$high_reads{$array[0]}="M2\t".$muts_included."\t".$remain_spikein."\t".$count_ref."\t".$muts_rep."\t".$line."\n";
								$reads_count_spk{$array[0]}=1;
							}    
						}
							
						#all M2 reads go through the slow mode check;
						$doubt_sam{$line."\n"}=1;
						$doubt_sam_ID{$array[0]}=$line."\n"; 
					}
				}
			}
		}
	}
}


#unpaired (high quality) reads were put at the end of the file.
foreach (keys %high_reads){ 
    print OFH_h $high_reads{$_};
}

#close the output file and complete.
close OFH_h;
close FH;

#write unpaired reads to fastq, 
# open(OFH_txt,'>'.$hq_header.".unpair.txt");
# open(OFH_fastq,'>'.$hq_header.".unpair.fastq");
# foreach (keys %high_reads){
#   my $line =$high_reads{$_};
#   my @array = split("\t",$line);
    
#   print OFH_txt $line;    
#   print OFH_fastq '@'.$array[5]."\n";
#   print OFH_fastq $array[14]."\n";
#   print OFH_fastq '+'."\n";
#   print OFH_fastq $array[15]."\n";    
# }
# close OFH_txt;
# close OFH_fastq;

## Generate paired SAM FILE. ####
foreach (keys %doubt_sam){
    # Only doubted reads are go through the slow mode validation.
    print OFH_SAM_SPK $_;
}
close OFH_SAM_SPK;


## SLOW MODE: use human reference genome on non-perfect match (match & fair) reads
my $slow_fq_file = $hq_header.".filter.1.fastq";
my $slow_fq_file_2 = $hq_header.".filter.2.fastq";
my %suspect_reads=();
my %spike_in_reads=();
if ($slow ne 'false'){
    #convert doubted sam to fastq inputs
    system('java -jar Refs/picard.jar SamToFastq I='.$sam_file." FASTQ=".$slow_fq_file." SECOND_END_FASTQ=".$slow_fq_file_2.
            " VALIDATION_STRINGENCY=LENIENT 2>>tmp/samtofastq.log");
    if (-s $slow_fq_file > 0){
        #mapping doubted reads to human reference
        open(FH,'bwa mem -v 0 -B '.$BWA_B.' -O '.$BWA_O.' -t '.$threads.' '.$REF_homo.' '.
            $slow_fq_file.' '.$slow_fq_file_2.' 2>>tmp/bwa_run_2.log|');
        # open(OFH,'>'.$slow_output);
        
        while(<FH>){
            #ignore header information
            if ($_ !~ /^@/){
                my $line=$_;
                chomp($line);
                my @array=split("\t");
                if ($line =~ /MD:Z:([\dATCG]+)\t/){
                    #read MD string
                    my $sum_length_ref= length($array[9]);
                    my $md_string = $1;
                    my $mut_ref = scalar(split(/[ACGT]/,$md_string));
                    
                    if ($md_string !~ /[ATCG]0[ATCG]/){
                        #mapping to human reference should not include dinucleotides;
                        #measuring the mapping length - Ref
                        my $cigar_ref = $array[5];
                        my $left_clip=0;
                        my $right_clip=0;
                        $left_clip = $1 if ($cigar_ref =~ /(\d+)S\d+M/); # 
                        $right_clip = $1 if ($cigar_ref =~ /\d+M(\d+)S/); # 
                        my $match_ref = $sum_length_ref-$left_clip-$right_clip;
                        
                        #measuring the mapping length - SIC 
                        my @array_spk = split("\t",$doubt_sam_ID{$array[0]});
                        my $sum_length_spk = length($array_spk[9]); 
                        my $cigar_spk = $array_spk[5];
                        my $left_clip=0;
                        my $right_clip=0;
                        $left_clip = $1 if ($cigar_spk =~ /(\d+)S\d+M/); # 
                        $right_clip = $1 if ($cigar_spk =~ /\d+M(\d+)S/); # 
                        my $match_spk = $sum_length_spk-$left_clip-$right_clip;
                                
                        #read MD string from SIC reads.
                        my $md_spk = $1 if ($doubt_sam_ID{$array[0]} =~ /MD:Z:([\dATCG]+)\t/);
                        my $mut_spk = scalar(split(/[ACGT]/,$md_spk));
                            
                        #determining suspected reads    
                        if ($match_spk < $match_ref -10 || 
                            $md_spk =~ /[ATCG]0[ATCG]/ || 
                            $mut_spk > $mut_ref || 
                            ($mut_ref ==0 && $cigar_ref =~ /^\d+M$/) ){
                            
                            if (!exists $suspect_reads{$array[0]}){
                                #append suspected reads via read ID
                                $suspect_reads{$array[0]}=1;
                            }
                        }
                    }else{
                        # If contained dinucleotides, kick it out from the original/native reads pool.
                        if (exists $suspect_reads{$array[0]}){
							delete $suspect_reads{$array[0]};
						}
                    }
                }
            }
        }
        close FH;
    }
    ## CLEAN UP 
    system('rm -f '.$slow_fq_file);
    system('rm -f '.$slow_fq_file_2);   
}

# Make output file
open(OFH,'>'.$hq_file);
open(OFH_suspect,'>'.$hq_file_suspected);

# header/summary:
print OFH "Total Statistics:\n";
print OFH "1. Total Output Reads:\t".scalar(keys %reads_count_total)."\n";
print OFH "2. Total Mapped Reads:\t".scalar(keys %reads_count_bwa)."\n";
print OFH "3. Total Spike-in Reads (Before Filtering):\t".scalar(keys %reads_count_spk)."\n";
print OFH "4. M4 reads: (Spike-in >=4):\t".$spk_count{'>=4'}."\n";
print OFH "5. M3 reads: (Spike-in =3):\t".$spk_count{'3'}."\n";
print OFH "6. M2 reads: (Spike-in =2):\t".$spk_count{'2'}."\n";
print OFH "7. M1 reads: (Spike-in =1):\t".$spk_count{'1'}."\n";

# high-quality reads
open(FH,$hq_file_pref);
while(<FH>){
    my $line = $_;
    my @array=split("\t",$line);
    if (exists $suspect_reads{$array[5]}){
        if ($suspect_reads{$array[5]} == 1){
            #remove the suspected reads.
            print OFH_suspect $line;
        }else{
            print OFH $line;
            $spike_in_reads{$array[5]}= -1;
        }
    }else{
        print OFH $line;
        $spike_in_reads{$array[5]}=1;
    }
}
close FH;   
close OFH;
close OFH_suspect;

open(OFH_SAM_SPK,'>'.$sam_file);
foreach (keys %reliable_sam){
    # all reads will then output as the sam file, it is before slow mode validation!
    print OFH_SAM_SPK $_;
}
close OFH_SAM_SPK;

#output summary in the std-out
print "Total Statistics:\n";
print "1. Total Output Reads:\t".scalar(keys %reads_count_total)."\n";
print "2. Total Mapped Reads:\t".scalar(keys %reads_count_bwa)."\n";
print "3. Total Spike-in Reads (Before Filtering):\t".scalar(keys %reads_count_spk)."\n";

## generate fastq files
if ($type eq 'all' || $type eq 'gzip'){
    print "Done. Now generating output fastq/bam files ... \n";
    print "Warning: defined output files already exists; will overwrite! \n" if (-f $spike_fastq_1);
    
    # initialize
    my $spike_fastq_1_sam = $spike_fastq_1;
    my $origin_fastq_1_sam = $origin_fastq_1;
    my $suspect_fastq_1_sam = $suspect_fastq_1;
    
    my $signal = 0;
    my $count_pair_1= 0;
    my $read_line_i = 0;
    
    # separating fastq file 1
    if($SAMPLE_1){  
        print ('Extracting Reads ... (Pair - 1)'."\n");
        open(FH,'gzip -cd '.$SAMPLE_1.'|');
    }elsif($BAM){
        print ('Extracting Reads from BAM ... '."\n");
        if ($keepBam){
            open(FH,$samtools_path.' view '.$BAM.'|');
            
            $spike_fastq_1_sam =~ s/\.bam/\.sam/g;
            $origin_fastq_1_sam =~ s/\.bam/\.sam/g;
            $suspect_fastq_1_sam =~ s/\.bam/\.sam/g;

        }else{
            open(FH,$samtools_path.' fastq '.$BAM.'|');
        }
    }
    if ($keepBam){
        open(OFH_spike,'>'.$spike_fastq_1_sam);
        open(OFH_origin,'>'.$origin_fastq_1_sam);
        open(OFH_suspect,'>'.$suspect_fastq_1_sam);
        
        open(FH_header, $samtools_path.' view -H '.$BAM.'|');
        while(<FH_header>){
            # output original sam header;
            print OFH_spike $_;
            print OFH_origin $_;
            print OFH_suspect $_;
        }
        close FH_header;
    }else{
        open(OFH_spike,'>'.$spike_fastq_1);
        open(OFH_origin,'>'.$origin_fastq_1);
        open(OFH_suspect,'>'.$suspect_fastq_1);
    }
    

    while(<FH>){
        #read data from the original fastq/bam files.
        my $line = $_;
        if($keepBam){
            my @array = split("\t",$line);
            if (exists $spike_in_reads{$array[0]}){
                $count_pair_1++;
                $signal = 1;
            }elsif (exists $suspect_reads{$array[0]}){
                $signal = -1;
            }else{               
                $signal = 0;
            }
            print OFH_spike $line if ($signal ==1);
            print OFH_origin $line if ($signal ==0);
            print OFH_suspect $line if ($signal ==-1);
        }else{
            if ($read_line_i % 4 ==0){ # only test on the head line
                # Casava 1.8 format || Illumina reads ends with /1 or /2
                if (/^@(\S+)\s(\S+)/ || /^@(\S+)\/[12]$/){
                    my $array=$1;
                    if (exists $spike_in_reads{$array}){
                        $count_pair_1++;
                        $signal = 1;
                    }elsif (exists $suspect_reads{$array}){
                        $signal = -1;
                    }else{
                        $signal = 0;
                    }
                }elsif(/^@(.+)/){
                    my $array=$1;
                    chomp($array);
                    if (exists $spike_in_reads{$array}){
                        $count_pair_1++;
                        $signal = 1;
                    }elsif (exists $suspect_reads{$array}){
                        $signal = -1;
                    }else{
                        $signal = 0;
                    }
                }
            }
            $read_line_i = $read_line_i + 1 ;
            print OFH_spike $line if ($signal ==1);
            print OFH_origin $line if ($signal ==0);
            print OFH_suspect $line if ($signal ==-1);
        }
    }
    close FH;
    close OFH_spike;
    close OFH_origin;
    close OFH_suspect;
    
    
    # If paired-end reads;
    if($SAMPLE_2 ne 'empty'){
        # separating fastq file 2
        print ('Extracting Reads... (Pair - 2)'."\n");
        open(FH,'gzip -cd '.$SAMPLE_2.'|');
        open(OFH_spike,'>'.$spike_fastq_2)  ;
        open(OFH_origin,'>'.$origin_fastq_2);
        open(OFH_suspect,'>'.$suspect_fastq_2);
        
        $signal = 0;
        # my $count_pair_2=0;
        $read_line_i = 0;
        while(<FH>){
            my $line = $_;
            if ($read_line_i % 4 ==0){ # head line
                # Casava 1.8 format || Illumina reads ends with /1 or /2
                if (/^@(\S+)\s(\S+)/ || /^@(\S+)\/[12]$/){
                    my $array=$1;
                    if (exists $spike_in_reads{$array}){
                        # $count_pair_2++;
                        $signal = 1;
                    }elsif (exists $suspect_reads{$array}){
                        $signal = -1;
                    }else{
                        $signal = 0;
                    }
                }
            }
            $read_line_i = $read_line_i + 1 ;
            print OFH_spike $line if ($signal ==1);
            print OFH_origin $line if ($signal ==0);
            print OFH_suspect $line if ($signal ==-1);
        }
        close FH;
        close OFH_spike;
        close OFH_origin;
        close OFH_suspect;
    }
    
    if($keepBam){
        # Convert sam to bam file.
        system($samtools_path.' view -b '.$spike_fastq_1_sam.' >'.$spike_fastq_1);
        system($samtools_path.' view -b '.$origin_fastq_1_sam.' >'.$origin_fastq_1);
        # Remove original sam file.
        system('rm -f '.$spike_fastq_1_sam);
        system('rm -f '.$origin_fastq_1_sam);
    }
    
    # End of the processing.
    print ('Done! In total '.$count_pair_1." Spike-in reads have been separated.\n");
}

######   Generating Gzip files ######
if ($type eq 'gzip'){
    print ('Gzip output files.. (may take long time...)'."\n");
    system('gzip -1 '.$origin_fastq_1);
    system('gzip -1 '.$suspect_fastq_1);
    system('gzip '.$spike_fastq_1);
    if($SAMPLE_2){
        system('gzip -1 '.$origin_fastq_2);
        system('gzip -1 '.$suspect_fastq_2);
        system('gzip '.$spike_fastq_2);
    }
    print ('Program finished normally.'."\n");
}

###### evaluation on the separation file ######
if ($eval){
    #measure insert size in separation reads    
    system('java -jar Refs/picard.jar CollectInsertSizeMetrics I='.$sam_file.' O=insert_spikein.txt H=insert_spikein.pdf M=0.5');
}

###### CLEAN UP ######
if ($clean eq 'true'){
    system('rm -f tmp/bwa_run.log'); # remove BWA log file 
    system('rm -f tmp/bwa_run_2.log'); # remove BWA log file 
    system('rm -f tmp/bwa_run_3.log'); # remove BWA log file 
    system('rm -f tmp/samtofastq.log'); # remove BWA log file 
}

## Main Program End

## Subroutines for Test mode                 
sub get_test_data {
    my ($test_mode) = @_;
    my $SAMPLE_1 = '';
    my $SAMPLE_2 = 'empty';
    my $BAM = '';
    my $hq_header = '';
    print("In Test Mode\n");
    
    if ($test_mode eq 'test_pos') {
        # test Spike-in Sample
        $SAMPLE_1='data/test_pos_2_1.fastq.gz';
        $SAMPLE_2='data/test_pos_2_2.fastq.gz';
        
        $hq_header = $SAMPLE_1; 
        $hq_header =~ s/_1\.fastq\.gz//g;
        $hq_header =~ s/.*\///g;
    }
    if ($test_mode eq 'test_neg') {     
        # test 'Negative' Sample
        $SAMPLE_1='data/test_neg_2_1.fastq.gz'; 
        $SAMPLE_2='data/test_neg_2_2.fastq.gz';

        $hq_header = $SAMPLE_1; 
        $hq_header =~ s/_1\.fastq\.gz//g;
        $hq_header =~ s/.*\///g;
    }
    if ($test_mode eq 'test_bam') {     
        # test 'Positive' Sample in BAM format
        $BAM='data/test_pos.bam'; 

        $hq_header = $BAM;
        $hq_header =~ s/\.bam/_bam/g;   
        $hq_header =~ s/.*\///g;
    }
    return ($SAMPLE_1, $SAMPLE_2, $BAM, $hq_header);
}

sub prepare_ref {
    my ($ref_prefix) = @_;
    print ("customized spike-in reference are used. (".$ref_prefix.") Check required files ... \n");
    # Create BWA INDEX if not existed.
    if (! -f $ref_prefix.'.fasta.bwt'){
        print ("BWA index files are not detected, try generating them automatically... \n");
        system('bwa index '.$ref_prefix.'.fasta '); 
        if (-f $ref_prefix.'.fasta.bwt'){
            print ("Successfully generate BWA index. \n");
        }else{
            print ("Failed, please check BWA availability ... \n");
            exit 0;
        }
    }else{
        print ("BWA INDEX found!\n");
    }
    # Create the $mut_ref file of new spike-in reference if not existed.
    $mut_REF=$ref_prefix.'.vcf';
    if (! -f $mut_REF){
        print ("No Mutation position file detected, try generating them automatically... \n");
        open(FH,$ref_prefix.'.fasta');
        open(OFH,'>', $mut_REF);
        # Initialize
        my $curr_gene_name = '';
        my $i = 0;
        
        while( <FH> ){
            if (/^>(\S+)/){
                # reset position to 0 when a new gene/chrom started.
                $curr_gene_name = $1;
                $i = 0;
            }else {
                my $string = $_;
                chomp($string);
                foreach (split //, $string){
                    $i = $i + 1;
                    if ($_ =~ /[ACGT]/){
                        print OFH $curr_gene_name."\t".$i."\n";
                    }
                }
            }
        }
        
        close FH;
        close OFH;
        
        if (-f $mut_REF){
            print("Success.\n");
        }else{
            print("Failed... please contact the developer.\n");
            exit 0;
        }
    }else{
        print ("Mutation position file found!\n");
    }
    return ($mut_REF);
}
