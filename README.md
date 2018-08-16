# SIST: Separate Accugenomics Spike-In reads
## Version 1.30

08/16/2018: optimized for UMI reads separation; 

## Introduction

Spike-In Separation Toolbox (SIST) is developed for separate Accugenomics sequencing spike-in controls (reads） from original sample reads. The separation is majorly based on the unique characteristics of spike-in controls, which contains several di-nucleotides inside the target region.  

This toolbox has been specificially designed and used in *Sequencing Quality Control 2* project - Working Group II for the purpose of target sequencing quality control. To read more details about SEQC2 project, please visit [here](https://www.fda.gov/ScienceResearch/BioinformaticsTools/MicroarrayQualityControlProject/ucm507935.htm)   

## Tool Manual (Help) Page
use `perl SIST.pl -help` to get help information.

#### **Note**
During to privacy concern, currently the reference of Accugenomics Spike-in Control are not released. This tool **CANNOT** be ran successfully without the specific refernece genome. If you are project related person, please contact [Leihong](mailto:leihong.wu@fda.hhs.gov) for further assistance.  

##Usage:
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

## Example 

`perl SIST.pl -test_pos` (embeded spike-in sample; Roche Batch2 S8)  
	fastq_1: /SEQC2/upload/171003_SEQC2_Hyperplus_SIC_170bp/sample4rep2_S8_L001_R1_001.fastq.gz  
	fastq_2: /SEQC2/upload/171003_SEQC2_Hyperplus_SIC_170bp/sample4rep2_S8_L001_R2_001.fastq.gz  

Equals to :  
`perl SIST.pl -fastq_1=data/test_pos_2_1.fastq.gz -fastq_2=data/test_pos_2_2.fastq.gz -O=test_pos_2 -type=match`
	
##	
	
`perl SIST.pl -test_neg` (embeded negative sample; Roche Batch2 S1)  
	fastq_1: /SEQC2/upload/171003_SEQC2_Hyperplus_SIC_170bp/sample1rep1_S1_L001_R1_001.fastq.gz  
	fastq_2: /SEQC2/upload/171003_SEQC2_Hyperplus_SIC_170bp/sample1rep1_S1_L001_R2_001.fastq.gz  

Equals to :  
`perl SIST.pl -fastq_1=data/test_neg_2_1.fastq.gz -fastq_2=data/test_neg_2_2.fastq.gz -O=test_neg_2 -type=match`
