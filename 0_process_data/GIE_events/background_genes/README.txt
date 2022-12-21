# This scripts extract mutation and copy number of information for the 100 randomly selected genes (see ../../background_genes/select_random_genes.ipynb) saved in the background_genes.tsv file

### Steps
# 1. Run the pipeline for Hartwig samples
conda activate snakemake
snakemake --profile slurm --snakefile pipeline_extract_genomic_info_backgronund_genes.py --config   o=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/hla_events/ i=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/whitelisted_samples.tsv d=hmf --drop-metadata   # only for whitelisted samples

conda activate global
python concat_reports_hmf_background.py # this is the concatenated report that is used in the manuscript analysis, a table with all whitelisted samples

# 2. Run the pipeline for PCAWG samples
conda activate snakemake
snakemake --profile slurm --snakefile pipeline_extract_genomic_info_backgronund_genes.py --config   o=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/hla_events/ i=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/whitelisted_samples.txt d=hmf --drop-metadata   # only for whitelisted samples


conda activate global
python concat_reports_pcawg_background.py # this is the concatenated report that is used in the manuscript analysis, a table with all whitelisted samples