#!/bin/bash
input_bam=$1
outpath=$2
file_name=$3

module load sambamcram/samtools/1.7
# create index
samtools index $1

# First need to extract the coordinates of the HLA locus
echo "samtools view $1 "6:29854528-32726735" -b > $outpath/${file_name}_HLA.bam"
samtools view $1 "6:29854528-32726735" -b > $outpath/${file_name}_HLA.bam

# create a bam with the chr format
echo ""
echo "samtools view -h  $outpath/${file_name}_HLA.bam| sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $outpath/${file_name}_HLA.tmp.bam"
samtools view -h  $outpath/${file_name}_HLA.bam| sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $outpath/${file_name}_HLA.tmp.bam
echo "rm $outpath/${file_name}_HLA.bam"
rm $outpath/${file_name}_HLA.bam

# Convert to hg38 coordinates, change the path to the location of the liftover binary
echo ""
echo "/hpc/local/CentOS7/cog/software/miniconda3/envs/global/bin/CrossMap.py bam /hpc/local/CentOS7/cog/software/xhla/hg19ToHg38.over.chain.gz $outpath/${file_name}_HLA.tmp.bam $outpath/${file_name}_HLA_hg38"
/hpc/local/CentOS7/cog/software/miniconda3/envs/global/bin/CrossMap.py bam /hpc/local/CentOS7/cog/software/xhla/hg19ToHg38.over.chain.gz $outpath/${file_name}_HLA.tmp.bam $outpath/${file_name}_HLA_hg38
echo "rm $outpath/${file_name}_HLA.tmp.bam"
rm $outpath/${file_name}_HLA.tmp.bam

# Run the HLA typing
echo ""
echo "singularity run -B $outpath/:/output/ /hpc/local/CentOS7/cog/software/xhla/xhla.sif  --sample_id ${file_name} --input_bam_path /output/${file_name}_HLA_hg38.sorted.bam  --output_path /output/"
mkdir -p $outpath/tmp/ # change the path the the location of your singularity container
singularity run -B $outpath/:/output/  -B $outpath/tmp/:/tmp/ --pwd /tmp/ /hpc/local/CentOS7/cog/software/xhla/xhla.sif  --sample_id ${file_name} --input_bam_path /output/${file_name}_HLA_hg38.sorted.bam  --output_path /output/
echo "rm $outpath/${file_name}_HLA_hg38.sorted.bam"
rm $outpath/${file_name}_HLA_hg38.sorted.bam
echo "rm $outpath/${file_name}_HLA_hg38.sorted.bai"
rm $outpath/${file_name}_HLA_hg38.sorted.bai
rm -rf $outpath/tmp



