# Scripts to process purple .vcf mutations into a tab-delimited format that can be used by dndscv and other tools
# Please note that this can only be reproduced when having access to the .vcf files of purple

output_path_hmf="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/"
output_path_pcawg="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/"


conda activate snakemake 

## 1. Run processing of mutations in Hartwig and PCAWG samples

snakemake --profile slurm  --snakefile run_vaf_extraction_hmf.py --config   o=$output_path_hmf/lilac/ i=$output_path_hmf/list_samples_hla.tsv  --drop-metadata --latency-wait 10

snakemake --profile slurm  --snakefile run_vaf_extraction_pcawg.py --config   o=$output_path_pcawg/lilac/ i=$output_path_pcawg/list_donors_hla_and_muts.txt  --drop-metadata --latency-wait 10
