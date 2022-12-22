# This scripts extract mutation and copy number of information for the 100 randomly selected genes (see ../../background_genes/select_random_genes.ipynb) saved in the background_genes.tsv file

output_path_hmf="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/" # change to your location
output_path_pcawg="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/" # change to your location



### Steps
# 1. Run the pipeline for Hartwig samples
conda activate snakemake
snakemake --profile slurm --snakefile pipeline_extract_genomic_info_backgronund_genes.py --config   o=$output_path_hmf/hla_events/ i=$output_path_hmf/whitelisted_samples.tsv d=hmf --drop-metadata   # only for whitelisted samples

conda activate global
python concat_reports_hmf_background.py # this is the concatenated report that is used in the manuscript analysis, a table with all whitelisted samples

# 2. Run the pipeline for PCAWG samples
conda activate snakemake
snakemake --profile slurm --snakefile pipeline_extract_genomic_info_backgronund_genes.py --config   o=$output_path_pcawg/hla_events/ i=$output_path_pcawg/whitelisted_samples.txt d=pcawg --drop-metadata   # only for whitelisted samples


conda activate global
python concat_reports_pcawg_background.py # this is the concatenated report that is used in the manuscript analysis, a table with all whitelisted samples