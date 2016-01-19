# Repeat filtering program

## parse_repeats.pl 

Perl script to run the repeat filtering algorithm starting from BAM files. The script performs the following steps:
1. Modifies the reads IDs of the right reads
2. Merges the right and left reads
3. Sorts the merged reads by ID
4. Filters alignments using filter_repeats_bam.py 

### Program usage

`
perl parse_repeats.pl -1 left.bam -2 right.bam [-o outpre -d distance -m/--multi -s/--nosingle -p cpus –n/--debug -S/--samold]
`

| -1 FILE		|	aligned left reads. BAM format. |
| -2 FILE		|	aligned right reads. BAM format.|
| -o TEXT		|	prefix for output files. Can include directory, e.g. ‘–o test/parsed’ will write files in directory ‘test’ 
				with prefix ‘parsed’. In its absence files will be output in current directory. |
| -d INT        	| maximum acceptable distance between mates. Includes the length of forward and reverse reads. Default=1000.|
| -m/--multi    	| allow multi-mapping reads. |
| -s/--nosingle		| skip singleton alignments. |
| -p INT        	| number of CPUs for merging and sorting steps. Default=1. |
| -S/--samold   	| use this option if samtools version < 0.1.19. |

## filter_repeats_bam.py

Python program to perform filtering. Briefly,
a. Calculates the distance between right and left reads and filters alignments for which distance < user-defined value. If a read pair maps too far, both reads are discarded.
b. If multiple alignments are acceptable, filters alignments with the minimum total # of mismatches for right+left read.
c. If only left or right read of a pair is mapped (singleton), the alignment with minimum # of mismatches are kept.
d. If multi-mapping reads are not allowed, and multiple alignments are obtained after filtering, the read pair or singleton read is discarded.

### Program usage
`
python filter_repeats_bam.py [-o <outfile>] [-d/--distance maxDist] [-n/--debug] [-m/--multi] [-s/--nosingle]
`

| -o FILE          |   output file. SAM format. By default, output goes to STDOUT. |
| -d/--distance INT  | maximum acceptable distance between mates. Includes the length of forward and reverse reads. Default=1000. |
| -n/--debug         | run in debug mode. Outputs 100 valid alignments to outfile and reports details about processed and parsed alignments. |
| -m/--multi         | allow multi-mapping reads |
| -s/--nosingle	   |	skip singleton alignments |

## Notes:
1. This script can be run by itself if the first three steps - modifying right read IDs, merging left and right reads, sort by read ID - have been performed already. 
2. It is designed to work with a stream so samtools output can be piped directly to this program.
3. Since the output is in SAM format, it can be compressed to BAM with samtools.

