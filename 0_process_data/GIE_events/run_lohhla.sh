#!/bin/bash

# this script is part of internal controls

input_bam_normal=$1
input_bam_tumor=$2
output_final_path=$3
sample_id=$4
purity_info=$5
hlas_path=$6

hlas_fasta="/hpc/local/CentOS7/cog/software/polysolver/data/abc_complete.fasta"
hlas_exon="/hpc/local/CentOS7/cog/software/polysolver/data/hla.dat"

output_path=$output_final_path/tmp

echo "output path in ... $output_path"

if [ ! -d $output_path ]; then

    mkdir -p $output_path

fi

path_bams=$output_path/bams_$sample_id
normal_bam=$path_bams/normal.bam
tumor_bam=$path_bams/tumor_$sample_id.bam

# first convert bam

if [ ! -d $path_bams ]; then

mkdir -p $path_bams

fi

module load sambamcram/samtools/1.7


# normal

#samtools view -h  $input_bam_normal | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $normal_bam # only for chr6
#echo "samtools view -h  $input_bam_normal | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $normal_bam"
cp $input_bam_normal $normal_bam
samtools index $normal_bam



# tumor
#samtools view -h  $input_bam_tumor | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $tumor_bam # only for chr6
#echo "samtools view -h  $input_bam_tumor | sed -e '/^@SQ/s/SN\:/SN\:chr/' -e '/^[^@]/s/\t/\tchr/2' | samtools view -bS - > $tumor_bam"
cp $input_bam_tumor $tumor_bam
samtools index $tumor_bam


# run lohhla
source  $(conda info --base)/etc/profile.d/conda.sh
conda activate LOHHLA
path_lohhla="/hpc/local/CentOS7/cog/software/lohhla/LOHHLAscript.R"


echo "Rscript $path_lohhla --patientId $sample_id --outputDir $output_path --normalBAMfile $normal_bam --BAMDir $path_bams --hlaPath $hlas_path --HLAfastaLoc $hlas_fasta --CopyNumLoc $purity_info --mappingStep TRUE --minCoverageFilter 5 --fishingStep TRUE --cleanUp TRUE --HLAexonLoc $hlas_exon --numMisMatch 2"
Rscript $path_lohhla --patientId $sample_id --outputDir $output_path --normalBAMfile $normal_bam --BAMDir $path_bams --hlaPath $hlas_path --HLAfastaLoc $hlas_fasta --CopyNumLoc $purity_info --mappingStep TRUE --minCoverageFilter 5 --fishingStep TRUE --cleanUp TRUE --HLAexonLoc $hlas_exon --numMisMatch 2

# clean the output

touch $output_final_path/${sample_id}.5.DNA.HLAlossPrediction_CI.xls # create to avoid error, will be overwritten
mv $output_path/Figures/*.pdf  $output_final_path
mv $output_path/*.xls  $output_final_path
rm -rf $output_path
#conda deactivate
