#!usr/bin/perl -w

# control script to run the repeat filtering algorithm

# 1. add /2 to right read IDs
# 2. merge with left reads
# 3. sort by read id
# 4. filter

# dependencies: samtools, python, perl
# filter reads from repeat regions before using this script

use strict;
use warnings;

my $currdir = $0;
$currdir =~ s/\/[^\/]+$//;

my $left = "";
my $right = "";
my $outpre = "";
my $dist = 1000;
my $multi = 0;
my $cpus = 1;
my $filter = $currdir."/filter_repeats_bam.py";
my $samold = 0;
foreach my $i (0..scalar(@ARGV)-1){
	if($ARGV[$i] eq "-1"){
		$left = $ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-2"){
		$right = $ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-o"){
		$outpre = $ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-d"){
		$dist = $ARGV[$i+1];
	}	
	elsif($ARGV[$i] eq "-m" || $ARGV[$i] eq "--multi"){
		$multi = 1;
	}
	elsif($ARGV[$i] eq "-p"){
		$cpus = $ARGV[$i+1];
	}
	elsif($ARGV[$i] eq "-s" || $ARGV[$i] eq "--samold"){
		$samold = 1;
	}
}

if($left eq "" || $right eq ""){
	print "Program usage:\n\tperl parse_repeats.pl -1 left -2 right"; 
	print "[-o outfile -d distance -m/--multi -p cpus -s/--samold]\n\n";
	die "Warning: Do not use -p option if samtools version < 0.1.19\n";
}

# add /2 to read IDs for right reads
my $right_id = "";
if($right =~ /(.+)\.(.+)$/){
	$right_id = $1."-2.".$2;
}
else {
	die "Right read file format not recognized!\n";
}

# changing read IDs of right reads
print "\n1. Changing read IDs of right reads\n";
my $cmd = "samtools view -h ".$right." | awk 'BEGIN{OFS=\"\\t\"}{if(\$0 !~ /^@/){\$1 = \$1 \"\/2\"} print \$0}' | samtools view -bS - > ".$right_id;
print $cmd,"\n";
my $system = `$cmd`;

# merging left and right reads
print "\n2. Merging left and right reads\n";
my $merged = "merged.bam";
if($outpre ne ""){
	$merged = $outpre."-".$merged;
}

if($samold){
	$cmd = "samtools merge ".$merged." ".$left." ".$right_id;
}
else{
	$cmd = "samtools merge -@ ".$cpus." ".$merged." ".$left." ".$right_id;
}
print $cmd,"\n";
$system = `$cmd`;

# sorting reads by ID
print "\n3. Sorting merged reads by ID\n";
my $sorted = "";
if($samold){
	$sorted = "merged-sorted";
}
else{
	$sorted = "merged-sorted.bam";
}

if($outpre ne ""){
	$sorted = $outpre."-".$sorted;
}

if($samold){
	$cmd = "samtools sort -n ".$merged." ".$sorted;
}
elsif($cpus > 1){
	$cmd = "samtools sort -\@ ".$cpus." -n -o ".$sorted." ".$merged;
}
else{
	$cmd = "samtools sort -n -o ".$sorted." ".$merged;
}
print $cmd,"\n";
$system = `$cmd`;

# run filtering program
print "\n4. Filtering alignments\n";
my $outfile = "merged-sorted-filtered.bam";
if($outpre ne ""){
	$outfile = $outpre."-".$outfile;
}

if(-e $filter){
	if($samold){
		$cmd = "samtools view -h ".$sorted.".bam | python ".$filter." -d ".$dist." | samtools view -bS - > ".$outfile;
	}
	else{
		$cmd = "samtools view -h ".$sorted." | python ".$filter." -d ".$dist." | samtools view -bS - > ".$outfile;
	}
}
else{
	die "Could not find ".$filter."! Stopping\n";
}
print $cmd,"\n";
$system = `$cmd`;
