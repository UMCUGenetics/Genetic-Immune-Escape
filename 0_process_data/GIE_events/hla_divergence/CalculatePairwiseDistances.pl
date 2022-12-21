#!/usr/bin/perl

use strict;
use warnings;

# Script written on 03.01.18 by Tobias Lenz
# Calculates pairwise distances between amino acid sequences, based on aa distance matrix
# The script expects an input file with aligned amino acid sequences (capital letters) of identical length in fasta format
# The script outputs the calculated distances in two different formats:
#    [fastaFile]_PairwiseDistanceMatrix.txt - A full-factorial matrix of all pairwise distances
#    [fastaFile]_PairwiseDistanceList.txt - A list of all pairwise distances

# Expected input from command line:
# Options:
# -d <string> | A full-factorial distance matrix for all amino acids 
#              (default: Grantham distance matrix following Grantham 1974 Science)
# -f <string> | Input file name with amino acid sequences in fasta format

my $inputLine = join(' ',@ARGV);
my $outputLine = "\n***********\nCommand line: " . $inputLine . "\n";

my $aaDistanceMatrix = "AAdistMatrix_Grantham.cnv";
if ($inputLine =~m/-d ([^ ]+)/) {
	$aaDistanceMatrix = $1;
}

my $fastaFile;
if ($inputLine =~m/-f ([^ ]+)/) {
	$fastaFile = $1;
} else {
	print "No input file name found!\n";
	print "\nExpected input from command line:\n";
	print "-d <string> | A full-factorial distance matrix for all amino acids (default: Grantham distance matrix following Grantham 1974 Science)\n";
	print "-f <string> | Input file name with amino acid sequences in fasta format\n";
	die "\nScript aborted!\n";
}
$outputLine = $outputLine . "Input sequence file: " . $fastaFile . "\nDistance matrix used: " . $aaDistanceMatrix . "\n***********\n";
print "$outputLine";
print "\nCalculating distances...\n";



# Read fasta file with aa sequences into hash
open (FASTA, '<', "$fastaFile") or die "Could not open input fasta file. $!";
my $ID;
my $sequence;
my $a;
my %SequenceArray = ();
while ($inputLine = <FASTA>) {
	chomp $inputLine;
	$inputLine =~ s/\r//g;
	if ($inputLine =~ m/^\s*$/) {next;}	#Ignore empty lines
	if ($inputLine =~ m/>/) {
		$sequence = "";
		($a, $ID) = split(/>/, $inputLine);
	}
	else {
		$sequence = $sequence . $inputLine;
		$SequenceArray{$ID} = $sequence;
	}
}
my @alleleList = sort (keys %SequenceArray);
close FASTA;

# Read aa distances into hash by combining all possible aa pairs
my $aa1;
my $aa2;
my $tempAA;
my $numberOfAAs;
my @aaList;
my $distance;
my %aaPairwiseDistance = ();
my $problemFound = 0;

open (AAMATRIX, '<', "$aaDistanceMatrix") or die "$!";
$inputLine = <AAMATRIX>; 		# read and process header line
chomp $inputLine;
while ($inputLine =~m/\t/) {
	($tempAA, $inputLine) = split(/\t/, $inputLine, 2);
	push(@aaList, $tempAA);
}
push(@aaList, $inputLine);

foreach $tempAA (@aaList) {
}
while ($inputLine = <AAMATRIX>) {
	chomp $inputLine;
	$a = 0;
	($aa1, $inputLine) = split(/\t/, $inputLine, 2);
	while ($inputLine =~m/\t/) {
		$a = $a + 1;
		($distance, $inputLine) = split(/\t/, $inputLine, 2);
		$aa2 = $aa1 . $aaList[$a];
		$aaPairwiseDistance{$aa2} = $distance;
	}
	$a = $a + 1;
	$aa2 = $aa1 . $aaList[$a];
	$aaPairwiseDistance{$aa2} = $inputLine;
}
close AAMATRIX;

# Calculate and output pairwise average distance between alleles
my $lengthAlleleArray = $#alleleList;
my $overallSequenceLength = length($SequenceArray{$alleleList[0]});
my $outputFile = $fastaFile . "_PairwiseDistanceMatrix.txt";
my $outputFilePairwise = $fastaFile . "_PairwiseDistanceList.txt";
my $sequenceLength;
my $distanceSum;
my $i;
my $j;
my $aaIndex;
my $problematicAA;

open (PAIRS, '>', "$outputFilePairwise") or die "$!\n";
print PAIRS "Allele_pairs\tDistance\n";			# print header line

open (MATRIX, '>', "$outputFile") or die "$!\n";
print MATRIX "Alleles";			# print header line
foreach $ID (@alleleList) {
	print MATRIX "\t$ID";
	if ($overallSequenceLength != length($SequenceArray{$ID})) {
		die "\n!!! Unequal sequence length detected. Script aborted!!!\n";
	}
}
print MATRIX "\n";

for ($i=0; $i <= $lengthAlleleArray; $i++) {
	print MATRIX "$alleleList[$i]";
	for ($j=0; $j <= $lengthAlleleArray; $j++) {
		$sequenceLength = $overallSequenceLength;
		$distanceSum = 0;
		if ($i ne $j) {				# Calculate sum of pairwise aa distance between ith and jth allele
			for ($aaIndex=0; $aaIndex < $sequenceLength; $aaIndex++) {
				$tempAA = substr($SequenceArray{$alleleList[$i]}, $aaIndex, 1) . substr($SequenceArray{$alleleList[$j]}, $aaIndex, 1);
				if ($tempAA =~ m/[^ARNDCQEGHILKMFPSTWYV]+/) {
					$problemFound = 1;
					$problematicAA = $tempAA;
					$sequenceLength = $sequenceLength - 1;
				} else {
					$distanceSum = $distanceSum + $aaPairwiseDistance{$tempAA};
				}
			}
		}
		$distance = $distanceSum / $sequenceLength;		# Average over sequence length
		print MATRIX "\t$distance";
		if ($j > $i) {
			print PAIRS "$alleleList[$i] _ $alleleList[$j]\t$distance\n";
		}
	}
	print MATRIX "\n";
}

close MATRIX;
close PAIRS;

print "\nCalculations successfully finished and output files created:\n";
print "\t\t$outputFile\n";
print "\t\t$outputFilePairwise\n\n";


if ($problemFound == 1) {
	print "\n\n!!! Warning: One or more non-standard characters encountered in sequence comparison. First occurrence: $problematicAA \n This might not be a problem, but you should be aware of what this means. Specific sites were ignored in affected pairwise comparisons!\n\n";
}
