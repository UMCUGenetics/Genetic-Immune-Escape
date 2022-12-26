# Genetic immune escape landscape in primary and metastatic cancer
Source code to reproduce the analysis of the paper. 
## :warning::warning:Repository under construction.


This repository contains all the data and analysis related to the preprint:https://www.biorxiv.org/content/10.1101/2022.02.23.481444v3

The Hartwig Medical Foundation tumor analytical pipeline can be found at: 

https://github.com/hartwigmedical/hmftools

LILAC can be found at: 

https://github.com/hartwigmedical/hmftools/tree/master/lilac


## Project Structure

This repository is structured as follows:

```shell

 ├──create_metadata_dataset.ipynb            # This is the script to create the metadata table (Supp. Table 3) and performing the blacklisting of samples
 
 ├─ 0_process_data/ # This folder contains the scripts necessary to process all raw files (.bam), somatic data (somatic mutations, CNVs) and expression (eg, T-cell infiltration) from both PCAWG and Hartwig samples
   └── mutations/ # To process .vcf from purple into a tab-delimited format
   └── GIE_events/ # To run LILAC (tumor-normal paired mode) and extract relevant patient specific information used in the manuscript
     └── background_genes/ # Scripts to extract the somatic alteration incidence of the 100 randomly selected genes (see background_genes/select_random_genes.ipynb)
   └── hla_typing_tools/ # Scripts to run LILAC, xHLA and polysolver independently for the germline and the tumor
├─ 1_benchmark # Scripts and notebooks for LILAC benchmark for HLA-I typing and LOH (Figure 1)
├─ 2_gie_prevalence # Scripts and notebooks for GIE prevalence in primary and metastatic tumors (Figure 2) and comparison (Figure 3)
├─ 3_positive_selection # Scripts and notebooks for HLA-I and non-HLA positive selection analyses (Figure 4 amd Figure 6)
   └── scripts_positive_selection/ # Scripts to run the positive selection analyses, both for mutations and CNVs across cancer types and pancancer in Hartwig and PCAWG
   └── analysis/ # Visualization notebooks
├── external_data/ # directory with external data used by the scripts
├── metadata/ # diretory with patient metadata used by the scripts, some files are not shared (or certain columns were removed) to avoid sharing patient-senstive data. 
├── README.md        # this file
└── environment_analysisA.yaml # Anaconda environment YAML file for specific analysis
```

## Data access

Please see https://www.biorxiv.org/content/10.1101/2022.06.17.496528v1 and
https://github.com/UMCUGenetics/primary-met-wgs-comparison for further description of the original data

### access PCAWG data
Somatic variant calls, gene driver lists, copy number profiles and other core data of the PCAWG cohort generated by the Hartwig analytical pipeline are available for download at https://dcc.icgc.org/releases/PCAWG/Hartwig. Researchers will need to apply to the ICGC data access compliance office (https://daco.icgc-argo.org) for the ICGC portion of the dataset. Similarly, users with authorized access can download the TCGA portion of the PCAWG dataset at https://icgc.bionimbus.org/files/5310a3ac-0344-458a-88ce-d55445540120. Additional information on accessing the data, including raw read files, can be found at https://docs.icgc.org/pcawg/data/.

### access Hartwig data
Metastatic WGS data and metadata from the Hartwig Medical Foundation are freely available for academic use through standardized procedures. Request forms can be found at https://www.hartwigmedicalfoundation.nl/en/data/data-acces-request/

### High-resolution HLA-I typing of 96 Hartwig patients

Raw sequencing data of the high-resolution HLA typing performed by GenDx can be downloaded via European Genome-phenome Archive (http://www.ebi.ac.uk/ega/) under accession number EGAD00001008643. 

### Supplementary Information

Supplementary Tables from the manuscript include relevant information to reproduce the analysis displayed in the manuscript. 

 
