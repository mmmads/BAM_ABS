Author: Saima Sultana Tithi and Hong Tran

####Introduction
BAM_ABS is a tool which simulates a Bayesian model that computes the posterior probability of mapping a multiread to each candidate genomic location, taking advantage of uniquely aligned reads. The inputs of this tool are a set of ambiguously mapped reads and a set of uniquely mapped reads from Bismark tool; and the output is the most probable genomic locations for those ambiguously mapped reads. The most probable genomic location is the location which has the highest calculated probability.

####System Requirement
This software package has been tested on Ubuntu 14.04 LTS. To run this program, user needs to have samtools, perl, bedtools, and g++ installed on his/her Linux/UNIX system. This program has been tested using g++ 4.8.4.

####Compilation
If you receive BAM_ABS as a compressed file, first decompress it. Then use the following commands to create the executable file:
```
make clean
make
```
For this command to work, the user needs g++ installed on his/her system. You can use the following command to install g++:
```
sudo apt-get install g++
```
	
####Execution
#####a) Pre-process the input data:
_Step 1_: Prior to execute this step, Samtools need to be installed on the system. After installing Samtools, run Samtools to get overlapped unique reads in sam format
* Input argument 1: Ambiguous reads in BED format (output of Bismark tool)
* Input argument 2: Unique reads in BAM format (output of Bismark tool)
* Output: Unique reads in SAM format with mapping quality greater than a given value
```
samtools view -L ambiguous_read_file.bed all_unique_reads.bam -q 20 > unique_reads.sam
```
The above command will only retain reads with MAQ(Mapping Quality) > 20 with no header

_Step 2_: Run the following command to get rid of duplicates from the unique reads
* Input: Unique reads in SAM format (output of step 1)
* Output: Unique reads with no duplicate in SAM format
```
sort -n -r -k3,3 -k4,4 -k5,5 unique_reads.sam|uniq -u > unique_reads_nodup.sam
```

_Step 3_: If Perl is not installed in the system, then prior to this step, Perl needs to be installed. Run Convert_to_bed.pl to convert unique read file to bed format.
* Input: Unique reads with no duplicate in SAM format (output of step 2)
* Output: Unique reads with no duplicate in BED format
```
perl Convert_to_bed.pl unique_reads_nodup.sam
```

_Setp 4_: Prior to execute this step, Bedtools need to be installed. After installing Bedtools, to get overlapped unique reads by using Bedtools, run the following command in the bedtools folder
* Input argument 1 (ambiguous_read_file.bed): Ambiguous reads in BED format
* Input argument 2 (unique_reads_nodup.bed): Unique reads with no duplicate in BED format (output of step 3)
* Output (unique_overlap_read_file.txt): All overlapping unique reads in txt format
```
./intersectBed -a ambiguous_read_file.bed -b unique_reads_nodup.bed -wb -wa > unique_overlap_read_file.txt
```

#####b) Score the multi-reads:
Run main.exe in BAM_ABS folder using the following command:
```
./main file.fa ambiguous_read_file unique_overlap_read_file.txt
```
Here,
* Input argument 1 (file.fa): The reference file in Fasta format
* Input argument 2 (ambiguous_read_file): The file containing all ambiguously mapped reads in Fastq format
* Input argument 3 (unique_overlap_read_file): The file containing all uniquely mapped reads which are overlapped with multi-reads or ambiguously mapped reads in txt format (output of step 4)
* Output (Reads_with_highest_probable_location.sam): Output file contains multi-reads along with the most probable genomic location in SAM format. This file only contains those multi-reads for which a probable genomic location can be calculated using our model.

####SNP and Methylation Rate
BAM_ABS folder also contains a file "snp_methylation_rate.txt". This file contains snp and methylation rate. If the user wants to change any rate, he or she needs to modify this file. Please do not delete this file or do not change the format of this file, only modify the rate part if necessary.

####Result
This tool will generate one output file: Reads_with_highest_probable_location.sam. Reads_with_highest_probable_location.sam contains multi-reads along with the most probable genomic location in SAM format.

####Example
* Input file:
 1. The reference file for mouse (in Fasta format): mm10.fa (you can download this file from http://hgdownload.cse.ucsc.edu/downloads.html#mouse)
 2. Multiread file (in Fastq format): L5_10_sample0.1_ambiguous_final
 3. Overlapping uniquely mapped reads: L5_sample0.1_10_unique_overlap.txt
* Output file: Multireads aligned to highest probable locations (in SAM format): Reads_with_highest_probable_location.sam

BAM_ABS command for the given example:
```
./main $BAM_ABS_Home$/mm10.fa $BAM_ABS_Home$/Sample_input_output/input/L5_10_sample0.1_ambiguous_final $BAM_ABS_Home$/Sample_input_output/input/L5_sample0.1_10_unique_overlap.txt
```
