# Scripts run the positive selection anlayses both for point mutations and CNVs (amplifications, deletions and LOH for HLA-I)
# Please note that this can only be reproduced when having access to the somatic files from purple (both for Hartwig and PCAWG)

output_path_hmf="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/hmf/"
output_path_pcawg="/hpc/cuppen/projects/P0020_genetics_immune_escape/large_scale_primary_met/processed/pcawg/"


conda activate snakemake 

## 1. Run dNdScv based somatic mutation positive selection analyses

snakemake --profile slurm --snakefile mutations/run_driver_pipeline_dndscv.py --config i=$output_path_hmf/hmf_by_ttype.json o=$output_path_hmf/positive_selection/ --drop-metadata --latency-wait 20 

snakemake --profile slurm --snakefile mutations/run_driver_pipeline_dndscv.py --config i=$output_path_pcawg/pcawg_by_ttype.json o=$output_path_pcawg/positive_selection/ --drop-metadata --latency-wait 20 


## 2.1 Run CNV positive selection analyses


snakemake --profile slurm --snakefile cnv/pipeline_pos_selection_cnv.py --config i=$output_path_hmf/hmf_by_ttype.json o=$output_path_hmf/positive_selection/cnv/ d=HMF --drop-metadata --latency-wait 20 

snakemake --profile slurm --snakefile cnv/pipeline_pos_selection_cnv.py --config i=$output_path_pcawg/pcawg_by_ttype.json o=$output_path_pcawg/positive_selection/cnv/ d=PCAWG--drop-metadata --latency-wait 20 


## 2.2. Parse output positive CNV selection analysis, the output will be a subfloder out of the output_directory (eg, positive_selection/cnv/processed/) 

conda activate global 

python cnv/process_output_LOH.py  --chunks_representation ../../external_data/regions_kb.tsv.gz --input_directory $output_path_hmf/positive_selection/cnv/
python cnv/process_output_LOH.py  --chunks_representation ../../external_data/regions_kb.tsv.gz --input_directory $output_path_pcawg/positive_selection/cnv/

python cnv/process_output_AMP.py  --chunks_representation ../../external_data/regions_kb.tsv.gz --input_directory $output_path_hmf/positive_selection/cnv/
python cnv/process_output_AMP.py  --chunks_representation ../../external_data/regions_kb.tsv.gz --input_directory $output_path_pcawg/positive_selection/cnv/

python cnv/process_output_DEL.py  --chunks_representation ../../external_data/regions_kb.tsv.gz --input_directory $output_path_hmf/positive_selection/cnv/
python cnv/process_output_DEL.py  --chunks_representation ../../external_data/regions_kb.tsv.gz --input_directory $output_path_pcawg/positive_selection/cnv/