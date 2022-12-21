## General remark
This package includes two perl scripts, two amino acid distance matrixes (one based on scoring scheme proposed by Grantham 1974 Science [default], and one uniform matrix for simple p-distance calculation), and two example input files. It also contains alignments of human MHC (HLA) alleles that can be used directly when working with human MHC/HLA (only common and well-documented HLA alleles following the CWD v2 catalog from Mack et al. 2013). The scripts and package are provided and updated by Tobias Lenz (lenz@post.harvard.edu), and were initially published in Pierini & Lenz 2018 MolBiolEvol (DOI:10.1093/molbev/msy116). Please cite that paper if you use anything from this package.


## File descriptions

# Script 1: CalculatePairwiseDistances.pl
The perl script CalculatePairwiseDistances.pl calculates pairwise distances between amino acid (aa) sequences, based on a provided aa distance matrix
The script expects an input file with aligned amino acid sequences (capital letters) of identical length in fasta format
The script outputs the calculated distances in two different formats:
  [fastaFile]_PairwiseDistanceMatrix.txt - A full-factorial matrix of all pairwise distances
  [fastaFile]_PairwiseDistanceList.txt - A list of all pairwise distances

Expected input from command line:
Options:
-d <string> | A full-factorial distance matrix for all amino acids 
             (default: Grantham distance matrix following Grantham 1974 Science)
-f <string> | Input file name with amino acid sequences in fasta format

Example command line call:
perl CalculatePairwiseDistances.pl -d AAdistMatrix_Grantham.cnv -f HLA_ClassII_CWDonly.fas

# Script 2: CalculateIndividualDivergence.pl
The perl script CalculateIndividualDivergence.pl calculates the sequence distance between an individual's alleles (average and sum), based on aa distance matrix, independent of the actual number of alleles (allowing for copy-number variation). The distance measure is normalized by sequence length for each pairwise comparison (i.e. the sum of the distance scores between two alleles is divided by the total length of the alignment of the two alleles (excluding sites containing non-amino acid characters). Average distance is calculated among all provided alleles, so the user can decide whether to include identical alleles (e.g. homozygous locus) in the calculation (provide both allele copies) or not (provide only one allele copy). The zero distance of two identical alleles will bias the calculation of the mean in multi-locus genotypes. This might or might not be intended, so this point needs to be considered by the user. 
The script expects an input file with individual IDs and a list of alleles found in that individual. Furthermore, the script requires a fasta file with aligned amino acid sequences (equal lengths) of all alleles that occur in those individuals, and an amino acid distance matrix (matrixes for Grantham distance and p-distance are provided).
The script outputs the calculated distances as a tab-delimited text file:
   [InputFile]_IndividualDivergence.txt - A list of divergences per individual

Expected input from command line:
Options:
-d <string> | A full-factorial distance matrix for all amino acids 
             (default: Grantham distance matrix following Grantham 1974 Science)
-f <string> | Input file name with amino acid sequences in fasta format
-i <string> | Name of tab-delimited input file with individual genotypes in the following format (using the same allele naming scheme as in the fasta alignment):
ID1	DQB10201	DQB10201
ID2	DQB10203	DQB10202
ID5	DQB10202	DQB10201	DQB10302
ID6	DQB10203	DQB10201	DQB10303	DQB10301
ID8	DQB10305
ID9	DQB10301	DQB10309

Example command line call to calculate Grantham allele divergence:
perl CalculateIndividualDivergence.pl -d AAdistMatrix_Grantham.cnv -f HLA_ClassII_CWDonly.fas -i Test_IndividualGenotypes.txt

Example command line call to calculate divergence based on simple p-distance:
perl CalculateIndividualDivergence.pl -d AAdistMatrix_Uniform.cnv -f HLA_ClassII_CWDonly.fas -i Test_IndividualGenotypes.txt


## Test run (with human MHC alleles):
Provided are example files for calculation of individual MHC class II allele divergence (here using human HLA-DQB1 alleles as example). Humans do not show copy-number-variation at HLA-DQB1, but some other species do, so this is to be seen only as an example for the general calculation of allele divergence across one or more MHC loci (with or without CNV). To run the test run, just use the example comman line above. This should give an output at shown in the provided example output files.


## Calculating allele divergence for human MHC/HLA alleles:
The classical loci in humans are not duplicated, so that only one (homozygous) or two (heterozygous) alleles per locus are applicable. Allele divergence can be calculated per locus or also across loci (makes only sense within class I OR class II, not across classes). See also Chowell et al. 2019 Nat Medicine (https://www.nature.com/articles/s41591-019-0639-4) for more information on such cross-locus allele divergence. For calculation of human HLA allele divergence, fasta files of 'common and well-documented' alleles are provided for both class I (HLA-A, -B, -C) and class II (HLA-DRB1, -DQB1) genes (except for a few alleles with incomplete sequence), and can be used as they are (after confirming that allele names are matched between alignments and individual genotype file). Provided are alignments of the entire sequence of the antigen-binding domains (class I: exon 2+3, class II: exon2), as well as alignments with only the antigen binding residues (ABS/PBR sites, following Björkman et la. 1987 Nature and Brown et al. 1988 Nature, respectively), for calculating allele divergence only at these residues.
The script is first calculating all pairwise distances before it uses these to calculate individual allele divergence. This pairwise calculation scales of course with the number of alleles in the alignment. If calculating distances for a small number of individuals with a limited allele pool, it might make sense to subset the alignment files to the required alleles only, in order to speed up the calculation.


## Disclaimer: I have written and tested these scripts to the best of my knowledge. However, I'm only human. Ultimately it is the user's responsibility to check correct calculations and I can not be held accountable for possible mistakes in this code.