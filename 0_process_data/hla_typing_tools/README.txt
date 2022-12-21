# Scripts to run LILAC, xHLA and Polysolver in tumor-only and germline-only mode, only HLA-I typing. 
# Please note that this can only be reproduced when having access to the raw .bam/.cram files (necessary to perform HLA-I typing) 

conda activate snakemake 

## 1. Run LILAC HLA-I typing in Hartwig

snakemake --profile slurm  --snakefile pipeline_run_LILAC_HMF_validation.py --config   o=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/lilac/ i=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/list_samples_hla.tsv --drop-metadata --latency-wait 10

## 2. Run xHLA and polysolver HLA-I typing in Hartwig

snakemake --profile slurm --snakefile run_other_tools_tumor_normal_HMF.py --config   o=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/hla_events/ i=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/list_samples_hla.tsv --drop-metadata  -np

# Important!!
# Note xHLA processing script is run_hla_typing_xHLA.sh and uses the xhla.sif singulairty container
# Note Polysolver uses the Polysolver4 polysolver4_new.sif singulairty container

## 3. Run LILAC HLA-I typing in PCAWG

snakemake --profile slurm  --snakefile pipeline_run_LILAC_pcawg_validation.py --config   o=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/lilac/ i=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/list_donors_hla_and_muts.txt --drop-metadata --latency-wait 10

## 4. Run xHLA and polysolver HLA-I typing in Hartwig

snakemake --snakefile pipeline_hla_events_PCAWG.py --config o=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/hla_events/ i=/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/list_donors_hla_and_muts.txt --profile slurm  --latency-wait 20