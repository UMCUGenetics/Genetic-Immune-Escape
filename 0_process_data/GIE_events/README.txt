# This are the scripts used to pre-process the raw data of the Harwtig and PCAWG patients. 
# Please note that this can only be reproduced when having access to the raw .bam/.cram files (necessary to perform HLA-I typing, somatic mutations) as well as to the somatic files processed by the Hartwig Medical Foundation pipeline

output_path_hmf="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/" # change to your location
output_path_pcawg="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/" # change to your location

# Steps

conda activate snakemake 

## 1. Run tumor-normal paired LILAC for Hartwig and PCAWG, perform HLA-I typing, tumor HLA-I copy number estimation and somatic mutations 

snakemake --profile slurm  --snakefile pipeline_run_LILAC_HMF.py --config   o=$output_path_hmf/lilac/ i=$output_path_hmf/list_samples_hla.tsv  --drop-metadata --latency-wait 10

snakemake --profile slurm  --snakefile pipeline_run_LILAC_pcawg.py --config   o=$output_path_pcawg/lilac/ i=$output_path_pcawg/list_donors_hla_and_muts.txt  --drop-metadata --latency-wait 10

## 2. Extract relevant GIE events and HLA-typing Harwtig samples

snakemake --profile slurm --snakefile pipeline_hla_events_HMF.py --config   o=$output_path_hmf/hla_events/ i=$output_path_hmf/list_samples_hla.tsv --drop-metadata  --latency-wait 10

## 3. Extract relevant GIE events and HLA-typing  PCAWG samples

snakemake --profile slurm --snakefile pipeline_hla_events_PCAWG.py --config o=$output_path_pcawg/hla_events/ i=$output_path_pcawg/list_donors_hla_and_muts.txt --profile slurm  --latency-wait 10

## 4. Process output of individual samples and prepare full table, the output of this scripts provide the input used in the manuscript analyses 

conda activate global 
python concat_reports_hmf.py
python concat_reports_pcawg.py

## 5. Extact immune infiltration Hartwig

conda activate snakemake
snakemake --profile slurm --snakefile pipeline_infiltration_HMF.py --config   o=$output_path_hmf/infiltration/ i=$output_path_hmf/list_samples_hla.tsv --drop-metadata  --latency-wait 10

conda activate global
python concat_reports_infiltration_hmf.py # this outputs the file hmf_infiltration.tsv.gz that is used in the manuscript analyses

## 6. Extact immune infiltration PCAWG

# RNA data extracted from https://dcc.icgc.org/releases/PCAWG/transcriptome/transcript_expression/pcawg.rnaseq.transcript.expr.tpm.tsv.gz
# ENSEMBL mapping from Transcript to Hugo symbol https://www.ensembl.org/biomart/ --> gene_ensembl_hugo.tsv
# Metadata conversion from aliquout_id to donor_id -> supp_table_1 from  https://www.biorxiv.org/content/10.1101/2022.06.17.496528v1

conda activate global
python pipeline_infiltration_PCAWG.py --rna_data pcawg.rnaseq.transcript.expr.tpm.tsv.gz --ensembl_mapping  gene_ensembl_hugo.tsv --metadata_pcawg metadata_pcawg.tsv --output_file $output_path_pcawg/rna/infiltration_pcawg.tsv.gz # this outputs the infiltration_pcawg.tsv.gz file that is used in the manuscript analyses





