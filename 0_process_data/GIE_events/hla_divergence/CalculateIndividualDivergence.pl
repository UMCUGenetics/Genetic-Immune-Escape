#!/usr/bin/perl

use strict;
use warnings;

# Script written on 01.06.18 and updated on 17.05.19 by Tobias Lenz (lenz@post.harvard.edu)
# Calculates average sequence distance between an individual's alleles, based on aa distance matrix, independent of the actual number of alleles (allowing for copy-number variation)
# The distance measure is normalized by sequence length for each pairwise comparison (i.e. the sum of the distance scores between two alleles is divided by the total length of the alignment of the two alleles (excluding sites containing non-amino acid characters)
# Average distance is calculated among all provided alleles, so the user can decide whether to include identical alleles (e.g. homozygous locus) in the calculation (provide both allele copies) or not (provide only one allele copy). The zero distance of two identical alleles will bias the calculation of the mean in multi-locus genotypes. This might or might not be intended.
# The script expects an input file with individual IDs and a list of alleles found in that individual
# Furthermore, the script requires a fasta file with aligned amino acid sequences (equal lengths) of all alleles that occur in those individuals, and an amino acid distance matrix (default is provided)
# The script outputs the calculated distances as a tab-delimited text file:
#    [InputFile]_IndividualDivergence.txt - A list of average divergences per individual

# Expected input from command line:
# Options:
# -d <string> | A full-factorial distance matrix for all amino acids 
#              (default: Grantham distance matrix following Grantham 1974 Science is provided with this package)
# -f <string> | Input file name with amino acid sequences in fasta format
# -i <string> | Name of tab-delimited input file with individual genotypes in the following format (column names and number of alleles are flexible):
#                 IDs	Allele1	Allele2	Allele3	Allele4
#                 ID1	DQB10201	DQB10201
#                 ID2	DQB10203	DQB10202
#                 ID5	DQB10202	DQB10201	DQB10302
#                 ID6	DQB10203	DQB10201	DQB10303	DQB10301
#                 ID8	DQB10305
#                 ID9	DQB10301	DQB10309

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
	print "No Fasta file name found!\n";
	print "\nExpected input from command line:\n";
	print "-d <string> | A full-factorial distance matrix for all amino acids (default: Grantham distance matrix following Grantham 1974 Science)\n";
	print "-f <string> | Input file name with amino acid sequences in fasta format\n";
	print "-i <stirng> | Name of tab-delimited input file with individual genotypes\n";
	print "-s <string> | Sample id\n";
	print "-o <string> | Output path \n";
	die "\nScript aborted!\n";
}

my $genotypeFile;
if ($inputLine =~m/-i ([^ ]+)/) {
	$genotypeFile = $1;
} else {
	print "No genotype file name found!\n";
	print "\nExpected input from command line:\n";
	print "-d <string> | A full-factorial distance matrix for all amino acids (default: Grantham distance matrix following Grantham 1974 Science)\n";
	print "-f <string> | Input file name with amino acid sequences in fasta format\n";
	print "-i <stirng> | Name of tab-delimited input file with individual genotypes\n";
	print "-s <string> | Sample id\n";
	print "-o <string> | Output path \n";
	die "\nScript aborted!\n";
}

my $output_path;
if ($inputLine =~m/-o ([^ ]+)/) {
	$output_path = $1;
} else {
	print "No output path found!!\n";
	print "\nExpected input from command line:\n";
	print "-d <string> | A full-factorial distance matrix for all amino acids (default: Grantham distance matrix following Grantham 1974 Science)\n";
	print "-f <string> | Input file name with amino acid sequences in fasta format\n";
	print "-i <stirng> | Name of tab-delimited input file with individual genotypes\n";
	print "-s <string> | Sample id\n";
	print "-o <string> | Output path \n";
	die "\nScript aborted!\n";
}

my $sample_id;
if ($inputLine =~m/-s ([^ ]+)/) {
	$sample_id = $1;
} else {
	print "No sample id found!!\n";
	print "\nExpected input from command line:\n";
	print "-d <string> | A full-factorial distance matrix for all amino acids (default: Grantham distance matrix following Grantham 1974 Science)\n";
	print "-f <string> | Input file name with amino acid sequences in fasta format\n";
	print "-i <stirng> | Name of tab-delimited input file with individual genotypes\n";
	print "-s <string> | Sample id\n";
	print "-o <string> | Output path \n";
	die "\nScript aborted!\n";
}

$outputLine = $outputLine . "Genotype file: " . $genotypeFile . "\nFasta file: " . $fastaFile . "\nDistance matrix used: " . $aaDistanceMatrix . "\n***********\n";
print "$outputLine";
print "\nCalculating distances...\n";



# Read fasta file with aa sequences into hash
print "$fastaFile";
open (FASTA, '<', "$fastaFile") or die "Aborted! Fasta file not found! $!";
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

open (AAMATRIX, '<', "$aaDistanceMatrix") or die "Aborted! Distance matrix file not found! $!";
$inputLine = <AAMATRIX>; 		# read and process header line
chomp $inputLine;
$inputLine =~ s/\r//g;
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

# Calculate pairwise average distance between alleles and store in hash
my $lengthAlleleArray = $#alleleList;
my $overallSequenceLength = length($SequenceArray{$alleleList[0]});
my $sequenceLength;
my $distanceSum;
my $i;
my $j;
my $aaIndex;
my $problematicAA;
my %allelePairs = ();




foreach $ID (@alleleList) {
	if ($overallSequenceLength != length($SequenceArray{$ID})) {
		die "\n!!! Unequal sequence length in fasta file detected. Script aborted!!!\n";
	}
}

for ($i=0; $i <= $lengthAlleleArray; $i++) {
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
					# print "Problem: $alleleList[$i] $alleleList[$j] $tempAA";   # Debugging
				} else {
					$distanceSum = $distanceSum + $aaPairwiseDistance{$tempAA};
				}
			}
		}
		$distance = $distanceSum / $sequenceLength;		# Average over sequence length
		$allelePairs{$alleleList[$i]}{$alleleList[$j]} = $distance;
		$allelePairs{$alleleList[$j]}{$alleleList[$i]} = $distance;
	}
}

if ($problemFound == 1) {
	print "\n\n!!! Warning: One or more non-standard characters encountered in sequence comparison. First occurrence: $problematicAA\n This might not be a problem, but you should be aware of what this means. Specific sites were ignored in affected pairwise comparisons!\n\n";
}


## Reading genotype file, calculate average divergence and write output
my $alleleNumber;
my $uniqueAlleleNumber;
my $tempAllele;
my @alleles;
my %alleleHash;
my $count = 0;
my $pairs;
# Open output file and write header line
my $outputFile = $output_path . "/" . $sample_id . ".IndividualDivergence.txt";
print $outputFile;
print "EEE\n";
open (OUTPUT, '>', "$outputFile") or die "Aborted! Could not create output file! $!";
print OUTPUT "ID\tAlleleNumber\tDivergence_Average\tDivergence_Sum\tAlleleList\n";

# Open genotype file and process individuals
open (GENOTYPES, '<', "$genotypeFile") or die "Aborted! Input file with individual genotypes not found! $!";
$inputLine = <GENOTYPES>;  # skip header line

while ($inputLine = <GENOTYPES>) {
	$problemFound = 0;
	$pairs = 0;
	@alleles = ();
	%alleleHash = ();
	$distanceSum = 0;
	$alleleNumber = 0;
	chomp $inputLine;
	$inputLine =~ s/\r//g;
	($ID, $inputLine) = split(/\t/, $inputLine, 2);

	while ($inputLine =~m/\t/) {
		($tempAllele, $inputLine) = split(/\t/, $inputLine, 2);
		if ($tempAllele ne "") {
			push(@alleles, $tempAllele);
			$alleleNumber++;
		}
	}
	if ($inputLine ne "") {
		push(@alleles, $inputLine);
		$alleleNumber++;
	}
	
	foreach $tempAllele (@alleles) {
		$alleleHash{$tempAllele}++;
	}
	$uniqueAlleleNumber = keys %alleleHash;
	
	if ($alleleNumber == 1) {
		print OUTPUT "$ID\t1\t0\t0\t$alleles[0]\n";
	} elsif ($alleleNumber > 1) {
		for ($i=0; $i < ($alleleNumber-1); $i++) {
			for ($j=$i+1; $j < $alleleNumber; $j++) {
				if (exists $allelePairs{$alleles[$i]}{$alleles[$j]}) {
					$distanceSum = $distanceSum + $allelePairs{$alleles[$i]}{$alleles[$j]};
					$pairs++;
				} else {
					$problemFound = 2;
				}
			}
		}
		if ($problemFound == 2) {
			print "\n!! Warning: One or more alleles of individual $ID were not found in the fasta file. Distances could not be calculated and were assigned NA!\n\n";
			$distance = "NA";
			$distanceSum = "NA";
		} else {
			$distance = $distanceSum / $pairs;
		}
		print OUTPUT "$ID\t$uniqueAlleleNumber\t$distance\t$distanceSum";
		foreach $tempAllele (@alleles) {
			print OUTPUT "\t$tempAllele";
		}
		print OUTPUT "\n";
	}
	$count++;
}
close GENOTYPES;



print "\nCalculations finished ($count genotypes processed) and output file created:\n";
print "\t$outputFile\n\n";


